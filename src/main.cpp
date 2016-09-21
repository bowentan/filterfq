#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <command_options.hpp>
#include <fastq_filter.hpp>
#include <quality_system.hpp>
#include <iostream>
#include <fstream>
#include <iterator>

using namespace boost::program_options;
using namespace boost::iostreams;
using namespace boost::filesystem;
using namespace boost::posix_time;
using namespace boost::gregorian;
using namespace command_options;
using namespace fastq_filter;
using namespace quality_system;
using namespace std;

string log_title() {return "[filterfq | " + to_simple_string(second_clock::local_time()) + "] ";}

int main(int argc, char* argv[]) {
    try {
        ptime start_time = second_clock::local_time();
        // input variables
        bool does_check_quality_sys;
        const string quality_sys[5] = {"Sanger", "Solexa", "Illumina 1.3+", "Illumina 1.5+", "Illumina 1.8+"};
        int raw_quality_sys;
        bool prefer_specified_raw_quality_sys;
        float max_base_N_rate;
        float min_ave_quality;
        int min_base_quality;
        float max_low_quality_rate;
        vector<path> raw_fq;

        // output variables
        int clean_quality_sys;
        path out_dir;
        string out_basename;
        vector<path> clean_fq;
        vector<path> dropped_fq;
        
        
        options_description generic("Gerneric options");
        generic.add_options()
            ("help,h", "produce help message")
        ;
        
        options_description param("Input parameters & files", options_description::m_default_line_length * 2, options_description::m_default_line_length);
        param.add_options()
            ("checkQualitySystem,c", bool_switch(&does_check_quality_sys), "only check quality system of the fastq file")
            ("baseNrate,N", value<float>(&max_base_N_rate) -> default_value(0.05), "maximum rate of \'N\' base allowed along a read")
            ("averageQuality,Q", value<float>(&min_ave_quality) -> default_value(0), "minimum average quality allowed along a read")
            ("perBaseQuality,q", value<int>(&min_base_quality) -> default_value(5), "minimum quality per base allowed along a read")
            ("lowQualityRate,r", value<float>(&max_low_quality_rate) -> default_value(0.5), "maximum low quality rate along a read")
            ("rawQualitySystem,s", value<int>(&raw_quality_sys), "specify quality system of raw fastq\n  0: Sanger\n  1: Solexa\n  2: Illumina 1.3+\n  3: Illumina 1.5+\n  4: Illumina 1.8+")
            ("preferSpecifiedRawQualitySystem,p", bool_switch(&prefer_specified_raw_quality_sys), "indicate that user prefers the given quality system to process")
            ("rawFastq,f", value< vector<path> >(&raw_fq) -> required() -> multitoken(), "raw fastq file(s) that need cleaned, required")
        ;

        options_description output("Output parameters & files", options_description::m_default_line_length * 2, options_description::m_default_line_length);
        output.add_options()
            ("cleanQualitySystem,S", value<int>(&clean_quality_sys) -> default_value(4), "specify quality system of cleaned fastq, the same as rawQualitySystem")
            ("outDir,O", value<path>(&out_dir) -> default_value(current_path().string()), "specify output directory, not used if cleanFastq is specified")
            ("outBasename,o", value<string>(&out_basename), "specify the basename for output file(s), required if outDir is specified")
            ("cleanFastq,F", value< vector<path> >(&clean_fq) -> multitoken(), "cleaned fastq file name(s), not used if outDir or outBasename is specified")
            ("droppedFastq,D", value< vector<path> >(&dropped_fq) -> multitoken(), "fastq file(s) containing reads that are filtered out")
        ;
        
        options_description desc("useage: filterfq <option>");
        desc.add(generic).add(param).add(output);

        variables_map vm; 
        store(command_line_parser(argc, argv).options(desc).run(), vm);
        if (vm.count("help")) {
            cout << setprecision(2) << desc << "\n";
            return 0;
        }
        check_option_independency(3, vm, "cleanFastq", "droppedFastq", "outBasename");
        check_option_independency(3, vm, "cleanFastq", "droppedFastq", "outDir");
        check_option_dependency(4, vm, "rawFastq", "checkQualitySystem", "cleanFastq", "droppedFastq", "outBasename");
        check_option_dependency(1, vm, "cleanFastq", "droppedFastq");
        check_option_dependency(1, vm, "droppedFastq", "cleanFastq");
        notify(vm);    
        
        for (vector<path>::iterator p = raw_fq.begin(); p != raw_fq.end(); p++) {
            try {
                *p = canonical(*p);
            }
            catch (filesystem_error& e) {
                cerr << "error: No such file or directory: " << e.path1().string() << endl;
                return 1;
            }
        }

        if (does_check_quality_sys) {
            cout << log_title() << "Welcome to filterfq!" << endl
                << log_title() << "INFO -- Parameters: filterfq ";
            for (variables_map::iterator i = vm.begin(); i != vm.end(); i++) {
                const variable_value& v = i -> second;
                if (v.value().type() == typeid(bool) && v.as<bool>()) {
                    cout << "--" << i -> first << " ";
                }
                else if (v.value().type() == typeid(vector<path>)) {
                    string op = i -> first;
                    if (op == "rawFastq") {
                        cout << "--" << op << " ";
                        for (path p : raw_fq) {
                            cout << p.string() << " ";
                        }
                    }
                }
            }
            cout << endl;
            int* checked_quality_sys = check_quality_system(raw_fq[0]);
            cout << log_title() << "INFO -- After checking 100,000 random reads, min quality code is \'" 
                << (char)checked_quality_sys[0] 
                << "\' and max quality code is \'" 
                << (char)checked_quality_sys[1] 
                << "\', the quality system is probably " 
                << quality_sys[checked_quality_sys[2]] 
                << "." << endl;
            ptime end_time = second_clock::local_time();
            time_duration dt = end_time - start_time;
            cout << log_title() << "INFO -- Process finished successfully! "
                << dt.seconds() << " seconds elapsed. Thank you for using filterfq!" << endl;
            return 0;
        }
        
        try {
            out_dir = canonical(out_dir);
        }
        catch (filesystem_error& e) {
            cerr << "error: No such file or directory: " << e.path1().string() << endl;
            return 1;
        }

        if (raw_fq.size() != clean_fq.size()) {
            cerr << "error: unequal numbers of raw and clean fastq: "
                << raw_fq.size() << " raw fastq and "
                << clean_fq.size() << " clean fastq." << endl;
            return 1;
        }
        else if (raw_fq.size() != dropped_fq.size()) {
            cerr << "error: unequal numbers of raw and dropped fastq: "
                << raw_fq.size() << " raw fastq and "
                << dropped_fq.size() << " dropped fastq." << endl;
            return 1;
        }
        else if (clean_fq.size() != dropped_fq.size()) {
            cerr << "error: unequal numbers of raw and dropped fastq: "
                << clean_fq.size() << " clean fastq and "
                << dropped_fq.size() << " dropped fastq." << endl;
            return 1;
        }
        
        if (vm.count("outBasename")) {
            if (raw_fq.size() == 1) {
                clean_fq.push_back(out_dir / path(out_basename + ".clean.fastq.gz"));
                dropped_fq.push_back(out_dir / path(out_basename + ".dropped.fastq.gz"));
            }
            else if (raw_fq.size() == 2) {
                for (int i = 0; i < raw_fq.size(); i++) {
                    clean_fq.push_back(out_dir / path(out_basename + "_" + to_string(i + 1) + ".clean.fastq.gz"));
                    dropped_fq.push_back(out_dir / path(out_basename + "_" + to_string(i + 1) + ".dropped.fastq.gz"));
                }
            }
        }
        else {
            for (int i = 0; i < clean_fq.size(); i++) {
                try {
                    clean_fq[i] = canonical(clean_fq[i].parent_path()) / clean_fq[i].filename();
                    dropped_fq[i] = canonical(dropped_fq[i].parent_path()) / dropped_fq[i].filename();
                }
                catch (filesystem_error& e) {
                    cerr << "error: No such file or directory: " << e.path1().string() << endl;
                    return 1;
                }
            }
        }

        cout << log_title() << "Welcome to filterfq!" << endl
            << log_title() << "INFO -- Filtering parameters: ";
        for (variables_map::iterator i = vm.begin(); i != vm.end(); i++) {
            const variable_value& v = i -> second;
            // if (v.value().type() == typeid(string)) {
            //     cout << "--" << i -> first << " " << v.as<string>() << " ";
            // }
            if (v.value().type() == typeid(float)) {
                cout << i -> first << "=" << v.as<float>() << " ";
            }
            else if (v.value().type() == typeid(int)) {
                cout << i -> first << "=" << v.as<int>() << " ";
            }
            else if (v.value().type() == typeid(bool)) {
                cout << i -> first << "=" << boolalpha << v.as<bool>() << " ";
            }
            // else if (!v.defaulted() && v.value().type() == typeid(path)) {
            //     cout << "--" << i -> first << " " << out_dir.string() << " ";
            // }
            // else if (v.value().type() == typeid(vector<path>)) {
            //     string op = i -> first;
            //     cout << "--" << op << " ";
            //     if (op == "rawFastq") {
            //         for (path p : raw_fq) {
            //             cout << p.string() << " ";
            //         }
            //     }
            //     else if (op == "cleanFastq") {
            //         for (path p : clean_fq) {
            //             cout << p.string() << " ";
            //         }
            //     }
            //     else if (op == "droppedFastq") {
            //         for (path p : dropped_fq) {
            //             cout << p.string() << " ";
            //         }
            //     }
            // }
        }
        cout << endl;

        cout << log_title() << "INFO -- The cleaned fastq files will be writen to ";
        for (vector<path>::const_iterator p = clean_fq.begin(); p != clean_fq.end(); p++) {
            if ((p + 1) != clean_fq.end()) {
                cout << (*p).string() << ", ";
            }
            else {
                cout << (*p).string() << ".";
            }
        }
        cout << endl;
        cout << log_title() << "INFO -- The dropped fastq files will be writen to ";
        for (vector<path>::const_iterator p = dropped_fq.begin(); p != dropped_fq.end(); p++) {
            if ((p + 1) != dropped_fq.end()) {
                cout << (*p).string() << ", ";
            }
            else {
                cout << (*p).string() << ".";
            }
        }
        cout << endl;

        int* checked_quality_sys = check_quality_system(raw_fq[0]);
        cout << log_title() << "INFO -- After checking 100,000 random reads, min quality code is \'" 
            << (char)checked_quality_sys[0] 
            << "\' and max quality code is \'" 
            << (char)checked_quality_sys[1] 
            << "\', the quality system is probably " 
            << quality_sys[checked_quality_sys[2]] 
            << "." << endl;
        if (prefer_specified_raw_quality_sys) {
            cout << log_title() << "WARN -- User prefered specified quality system "
                << quality_sys[raw_quality_sys]
                << ". The program will treat quality codes correspondingly." << endl;
        }
        else {
            raw_quality_sys = checked_quality_sys[2];
            cout << log_title() << "INFO -- The program will treat quality codes in accordance to "
                << quality_sys[raw_quality_sys] << "." << endl;
        }
        
        if (raw_quality_sys == clean_quality_sys) {
            cout << log_title() << "INFO -- All quality codes will remain the same as they were in "
                << quality_sys[raw_quality_sys] << "." << endl;
        }
        else {
            cout << log_title() << "INFO -- All quality codes will be converted to the corresponding codes in "
                << quality_sys[clean_quality_sys] << "." << endl;
        }

        switch(raw_fq.size()) {
            case 1: 
                {
                    cout << log_title() << "INFO " 
                        << setw(12) << "Processed" << " | "
                        << setw(12) << "Filtered" << " | "
                        << setw(12) << "Ratio(%)" << " | "
                        << setw(12) << "High" << " | "
                        << setw(12) << "Ratio(%)" << " | "
                        << setw(12) << "Ratio(%)" << " | "
                        << setw(12) << "Low average" << " | "
                        << setw(12) << "Ratio(%)" << " | "
                        << setw(12) << "Ratio(%)" << " | "
                        << setw(12) << "High low-" << " | "
                        << setw(12) << "Ratio(%)" << " | "
                        << setw(12) << "Ratio(%)" << " | " 
                        << endl;
                    cout << log_title() << "INFO "
                        << setw(12) << "reads" << " | "
                        << setw(12) << "reads" << " | "
                        << setw(12) << "filtered" << " | "
                        << setw(12) << "N-rate" << " | "
                        << setw(12) << "in processed" << " | "
                        << setw(12) << "in filtered" << " | "
                        << setw(12) << "quality" << " | "
                        << setw(12) << "in processed" << " | "
                        << setw(12) << "in filtered" << " | "
                        << setw(12) << "quality-rate" << " | "
                        << setw(12) << "in processed" << " | "
                        << setw(12) << "in filtered" << " | "
                        << endl;
                    cout << log_title() << "INFO ";
                    for (int i = 0; i < 15 * 12; i++) {
                        cout << "-";
                    }
                    cout << endl;

                    ifstream infq(raw_fq[0].c_str(), ios_base::in | ios_base::binary);
                    filtering_istream infq_decompressor;
                    infq_decompressor.push(gzip_decompressor());
                    infq_decompressor.push(infq);

                    ofstream clean_outfq(clean_fq[0].c_str(), ios_base::out | ios_base::binary);
                    filtering_ostream clean_outfq_compressor;
                    clean_outfq_compressor.push(gzip_compressor());
                    clean_outfq_compressor.push(clean_outfq);

                    ofstream dropped_outfq(dropped_fq[0].c_str(), ios_base::out | ios_base::binary);
                    filtering_ostream dropped_outfq_compressor;
                    dropped_outfq_compressor.push(gzip_compressor());
                    dropped_outfq_compressor.push(dropped_outfq);
                    
                    int n_read_processed = 0,
                        n_read_high_N_rate = 0,
                        n_read_low_ave_quality = 0,
                        n_read_high_low_quality_rate = 0,
                        n_read_filtered = 0;
                    bool is_filtered;
                    string read_id_line, read_line, plus_line, quality_line;
                    while (getline(infq_decompressor, read_id_line)) {
                        getline(infq_decompressor, read_line);
                        getline(infq_decompressor, plus_line);
                        getline(infq_decompressor, quality_line);
                        n_read_processed++;
                        
                        is_filtered = false;
                        // cout << read_id_line << endl
                        //     << read_line << endl 
                        //     << plus_line << endl
                        //     << quality_line << endl;
                        // cout << get_base_N_rate(read_line) << " "
                        //     << get_average_quality(quality_line, raw_quality_sys) << " "
                        //     << get_low_quality_rate(quality_line, raw_quality_sys, base_quality) << endl;
                        if (get_base_N_rate(read_line) > max_base_N_rate) {
                            n_read_high_N_rate++;
                            if (!is_filtered) {
                                n_read_filtered++;
                                is_filtered = true;
                            }
                        }
                        if (get_average_quality(quality_line, raw_quality_sys) < min_ave_quality) {
                            n_read_low_ave_quality++;
                            if (!is_filtered) {
                                n_read_filtered++;
                                is_filtered = true;
                            }
                        }
                        if (get_low_quality_rate(quality_line, raw_quality_sys, min_base_quality) > max_low_quality_rate) {
                            n_read_high_low_quality_rate++;
                            if (!is_filtered) {
                                n_read_filtered++;
                                is_filtered = true;
                            }
                        }
                        
                        if (n_read_processed % 50000 == 0) {
                            cout << log_title() << "INFO " 
                                << fixed
                                << setw(12) << setprecision(6) << n_read_processed << " | "
                                << setw(12) << setprecision(6) << n_read_filtered << " | "
                                << setw(12) << setprecision(6) << n_read_filtered * 100.0 / n_read_processed << " | "
                                << setw(12) << setprecision(6) << n_read_high_N_rate << " | "
                                << setw(12) << setprecision(6) << n_read_high_N_rate * 100.0 / n_read_processed << " | "
                                << setw(12) << setprecision(6) << n_read_high_N_rate * 100.0 / n_read_filtered << " | "
                                << setw(12) << setprecision(6) << n_read_low_ave_quality << " | "
                                << setw(12) << setprecision(6) << n_read_low_ave_quality * 100.0 / n_read_processed << " | "
                                << setw(12) << setprecision(6) << n_read_low_ave_quality * 100.0 / n_read_filtered << " | "
                                << setw(12) << setprecision(6) << n_read_high_low_quality_rate << " | "
                                << setw(12) << setprecision(6) << n_read_high_low_quality_rate * 100.0 / n_read_processed << " | "
                                << setw(12) << setprecision(6) << n_read_high_low_quality_rate * 100.0 / n_read_filtered << " | "
                                << endl;
                        }

                        if (!is_filtered) {
                            if (raw_quality_sys != clean_quality_sys) {
                                quality_system_convert(quality_line, raw_quality_sys, clean_quality_sys);
                            }
                            // cout << get_base_N_rate(read_line) << " "
                            //     << get_average_quality(quality_line, clean_quality_sys) << " "
                            //     << get_low_quality_rate(quality_line, clean_quality_sys, base_quality) << endl;
                            clean_outfq_compressor << read_id_line << endl
                                << read_line << endl
                                << plus_line << endl
                                << quality_line << endl;
                        }
                        else {
                            dropped_outfq_compressor << read_id_line << endl
                                << read_line << endl
                                << plus_line << endl
                                << quality_line << endl;
                        }
                    }
                    close(infq_decompressor, ios_base::in);
                    close(clean_outfq_compressor, ios_base::out);
                    close(dropped_outfq_compressor, ios_base::out);
                    cout << log_title() << "INFO ";
                    for (int i = 0; i < 15 * 12; i++) {
                        cout << "-";
                    }
                    cout << endl;
                    cout << log_title() << "INFO " 
                        << fixed
                        << setw(12) << setprecision(6) << n_read_processed << " | "
                        << setw(12) << setprecision(6) << n_read_filtered << " | "
                        << setw(12) << setprecision(6) << n_read_filtered * 100.0 / n_read_processed << " | "
                        << setw(12) << setprecision(6) << n_read_high_N_rate << " | "
                        << setw(12) << setprecision(6) << n_read_high_N_rate * 100.0 / n_read_processed << " | "
                        << setw(12) << setprecision(6) << n_read_high_N_rate * 100.0 / n_read_filtered << " | "
                        << setw(12) << setprecision(6) << n_read_low_ave_quality << " | "
                        << setw(12) << setprecision(6) << n_read_low_ave_quality * 100.0 / n_read_processed << " | "
                        << setw(12) << setprecision(6) << n_read_low_ave_quality * 100.0 / n_read_filtered << " | "
                        << setw(12) << setprecision(6) << n_read_high_low_quality_rate << " | "
                        << setw(12) << setprecision(6) << n_read_high_low_quality_rate * 100.0 / n_read_processed << " | "
                        << setw(12) << setprecision(6) << n_read_high_low_quality_rate * 100.0 / n_read_filtered << " | "
                        << endl;
                    ptime end_time = second_clock::local_time();
                    time_duration dt = end_time - start_time;
                    cout << log_title() << "INFO -- Process finished successfully! "
                        << dt.seconds() << " seconds elapsed. Thank you for using filterfq!" << endl;
                    break;
                }
            case 2:
                {
                    cout << log_title() << "INFO " 
                        << setw(12) << "Processed" << " | "
                        << setw(12) << "Filtered" << " | "
                        << setw(12) << "Ratio(%)" << " | "
                        << setw(12) << "High" << " | "
                        << setw(12) << "Ratio(%)" << " | "
                        << setw(12) << "Ratio(%)" << " | "
                        << setw(12) << "Low average" << " | "
                        << setw(12) << "Ratio(%)" << " | "
                        << setw(12) << "Ratio(%)" << " | "
                        << setw(12) << "High low-" << " | "
                        << setw(12) << "Ratio(%)" << " | "
                        << setw(12) << "Ratio(%)" << " | " 
                        << endl;
                    cout << log_title() << "INFO "
                        << setw(12) << "pairs" << " | "
                        << setw(12) << "pairs" << " | "
                        << setw(12) << "filtered" << " | "
                        << setw(12) << "N-rate" << " | "
                        << setw(12) << "in processed" << " | "
                        << setw(12) << "in filtered" << " | "
                        << setw(12) << "quality" << " | "
                        << setw(12) << "in processed" << " | "
                        << setw(12) << "in filtered" << " | "
                        << setw(12) << "quality-rate" << " | "
                        << setw(12) << "in processed" << " | "
                        << setw(12) << "in filtered" << " | "
                        << endl;
                    cout << log_title() << "INFO ";
                    for (int i = 0; i < 15 * 12; i++) {
                        cout << "-";
                    }
                    cout << endl;

                    // iostream for the first fastq
                    ifstream infq1(raw_fq[0].c_str(), ios_base::in | ios_base::binary);
                    filtering_istream infq1_decompressor;
                    infq1_decompressor.push(gzip_decompressor());
                    infq1_decompressor.push(infq1);

                    ofstream clean_outfq1(clean_fq[0].c_str(), ios_base::out | ios_base::binary);
                    filtering_ostream clean_outfq1_compressor;
                    clean_outfq1_compressor.push(gzip_compressor());
                    clean_outfq1_compressor.push(clean_outfq1);

                    ofstream dropped_outfq1(dropped_fq[0].c_str(), ios_base::out | ios_base::binary);
                    filtering_ostream dropped_outfq1_compressor;
                    dropped_outfq1_compressor.push(gzip_compressor());
                    dropped_outfq1_compressor.push(dropped_outfq1);
                    
                    // iostream for the second fastq
                    ifstream infq2(raw_fq[1].c_str(), ios_base::in | ios_base::binary);
                    filtering_istream infq2_decompressor;
                    infq2_decompressor.push(gzip_decompressor());
                    infq2_decompressor.push(infq2);

                    ofstream clean_outfq2(clean_fq[1].c_str(), ios_base::out | ios_base::binary);
                    filtering_ostream clean_outfq2_compressor;
                    clean_outfq2_compressor.push(gzip_compressor());
                    clean_outfq2_compressor.push(clean_outfq2);

                    ofstream dropped_outfq2(dropped_fq[1].c_str(), ios_base::out | ios_base::binary);
                    filtering_ostream dropped_outfq2_compressor;
                    dropped_outfq2_compressor.push(gzip_compressor());
                    dropped_outfq2_compressor.push(dropped_outfq2);

                    int n_pair_processed = 0,
                        n_pair_high_N_rate = 0,
                        n_pair_low_ave_quality = 0,
                        n_pair_high_low_quality_rate = 0,
                        n_pair_filtered = 0;
                    bool is_filtered;
                    string read_id_line1, read_line1, plus_line1, quality_line1,
                           read_id_line2, read_line2, plus_line2, quality_line2;
                    while (getline(infq1_decompressor, read_id_line1) && getline(infq2_decompressor, read_id_line2)) {
                        getline(infq1_decompressor, read_line1);
                        getline(infq1_decompressor, plus_line1);
                        getline(infq1_decompressor, quality_line1);
                        getline(infq2_decompressor, read_line2);
                        getline(infq2_decompressor, plus_line2);
                        getline(infq2_decompressor, quality_line2);
                        n_pair_processed++;
                        
                        is_filtered = false;
                        
                        if (get_base_N_rate(read_line1) > max_base_N_rate || get_base_N_rate(read_line2) > max_base_N_rate) {
                            n_pair_high_N_rate++;
                            if (!is_filtered) {
                                n_pair_filtered++;
                                is_filtered = true;
                            }
                        }
                        if (get_average_quality(quality_line1, raw_quality_sys) < min_ave_quality || get_average_quality(quality_line2, raw_quality_sys) < min_ave_quality) {
                            n_pair_low_ave_quality++;
                            if (!is_filtered) {
                                n_pair_filtered++;
                                is_filtered = true;
                            }
                        }
                        if (get_low_quality_rate(quality_line1, raw_quality_sys, min_base_quality) > max_low_quality_rate || get_low_quality_rate(quality_line2, raw_quality_sys, min_base_quality) > max_low_quality_rate) {
                            n_pair_high_low_quality_rate++;
                            if (!is_filtered) {
                                n_pair_filtered++;
                                is_filtered = true;
                            }
                        }

                        if (n_pair_processed % 50000 == 0) {
                            cout << log_title() << "INFO " 
                                << fixed
                                << setw(12) << setprecision(6) << n_pair_processed << " | "
                                << setw(12) << setprecision(6) << n_pair_filtered << " | "
                                << setw(12) << setprecision(6) << n_pair_filtered * 100.0 / n_pair_processed << " | "
                                << setw(12) << setprecision(6) << n_pair_high_N_rate << " | "
                                << setw(12) << setprecision(6) << n_pair_high_N_rate * 100.0 / n_pair_processed << " | "
                                << setw(12) << setprecision(6) << n_pair_high_N_rate * 100.0 / n_pair_filtered << " | "
                                << setw(12) << setprecision(6) << n_pair_low_ave_quality << " | "
                                << setw(12) << setprecision(6) << n_pair_low_ave_quality * 100.0 / n_pair_processed << " | "
                                << setw(12) << setprecision(6) << n_pair_low_ave_quality * 100.0 / n_pair_filtered << " | "
                                << setw(12) << setprecision(6) << n_pair_high_low_quality_rate << " | "
                                << setw(12) << setprecision(6) << n_pair_high_low_quality_rate * 100.0 / n_pair_processed << " | "
                                << setw(12) << setprecision(6) << n_pair_high_low_quality_rate * 100.0 / n_pair_filtered << " | "
                                << endl;
                        }
                        
                        if (!is_filtered) {
                            if (raw_quality_sys != clean_quality_sys) {
                                quality_system_convert(quality_line1, raw_quality_sys, clean_quality_sys);
                                quality_system_convert(quality_line2, raw_quality_sys, clean_quality_sys);
                            }
                            clean_outfq1_compressor << read_id_line1 << endl
                                << read_line1 << endl
                                << plus_line1 << endl
                                << quality_line1 << endl;
                            clean_outfq2_compressor << read_id_line2 << endl
                                << read_line2 << endl
                                << plus_line2 << endl
                                << quality_line2 << endl;
                        }
                        else {
                            dropped_outfq1_compressor << read_id_line1 << endl
                                << read_line1 << endl
                                << plus_line1 << endl
                                << quality_line1 << endl;
                            dropped_outfq2_compressor << read_id_line2 << endl
                                << read_line2 << endl
                                << plus_line2 << endl
                                << quality_line2 << endl;
                        }
                    }
                    close(infq1_decompressor, ios_base::in);
                    close(infq2_decompressor, ios_base::in);
                    close(clean_outfq1_compressor, ios_base::out);
                    close(clean_outfq2_compressor, ios_base::out);
                    close(dropped_outfq1_compressor, ios_base::out);
                    close(dropped_outfq2_compressor, ios_base::out);
                    cout << log_title() << "INFO ";
                    for (int i = 0; i < 15 * 12; i++) {
                        cout << "-";
                    }
                    cout << endl;
                    cout << log_title() << "INFO " 
                        << fixed
                        << setw(12) << setprecision(6) << n_pair_processed << " | "
                        << setw(12) << setprecision(6) << n_pair_filtered << " | "
                        << setw(12) << setprecision(6) << n_pair_filtered * 100.0 / n_pair_processed << " | "
                        << setw(12) << setprecision(6) << n_pair_high_N_rate << " | "
                        << setw(12) << setprecision(6) << n_pair_high_N_rate * 100.0 / n_pair_processed << " | "
                        << setw(12) << setprecision(6) << n_pair_high_N_rate * 100.0 / n_pair_filtered << " | "
                        << setw(12) << setprecision(6) << n_pair_low_ave_quality << " | "
                        << setw(12) << setprecision(6) << n_pair_low_ave_quality * 100.0 / n_pair_processed << " | "
                        << setw(12) << setprecision(6) << n_pair_low_ave_quality * 100.0 / n_pair_filtered << " | "
                        << setw(12) << setprecision(6) << n_pair_high_low_quality_rate << " | "
                        << setw(12) << setprecision(6) << n_pair_high_low_quality_rate * 100.0 / n_pair_processed << " | "
                        << setw(12) << setprecision(6) << n_pair_high_low_quality_rate * 100.0 / n_pair_filtered << " | "
                        << endl;
                    ptime end_time = second_clock::local_time();
                    time_duration dt = end_time - start_time;
                    cout << log_title() << "INFO -- Process finished successfully! "
                        << dt.seconds() << " seconds elapsed. Thank you for using filterfq!" << endl;
                    break;
                }
        }
    }
    catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        return 1;
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
    }

    return 0;
}
