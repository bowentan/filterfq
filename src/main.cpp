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
        // input variables
        const string quality_sys[5] = {"Sanger", "Solexa", "Illumina 1.3+", "Illumina 1.5+", "Illumina 1.8+"};
        int raw_quality_sys;
        bool prefer_specified_raw_quality_sys;
        float base_N_rate;
        float ave_quality;
        int base_quality;
        float low_quality_rate;
        vector<path> raw_fq;

        // output variables
        int clean_quality_sys;
        path out_dir;
        string out_basename;
        vector<path> clean_fq;
        
        options_description desc("useage: filterfq");
        desc.add_options()
            // input
            ("help,h", "produce help message")
            ("baseNrate,N", value<float>(&base_N_rate) -> default_value(0.05), "maximum rate of \"N\" base allowed along a read")
            ("averageQuality,Q", value<float>(&ave_quality), "minimum average quality allowed along a read")
            ("perBaseQuality,q", value<int>(&base_quality) -> default_value(5), "minimum quality per base allowed along a read")
            ("lowQualityRate,r", value<float>(&low_quality_rate) -> default_value(0.5), "low quality rate along a read")
            ("rawQualitySystem,s", value<int>(&raw_quality_sys), "specify quality system of raw fastq")
            ("preferSpecifiedRawQualitySystem,p", bool_switch(&prefer_specified_raw_quality_sys), "indicate that user prefers the given quality system to process")
            ("rawFastq,f", value< vector<path> >(&raw_fq) -> multitoken(), "raw fastq file(s) that need cleaned, required")
            // output
            ("cleanQualitySystem,S", value<int>(&clean_quality_sys) -> default_value(4), "quality system of cleaned fastq")
            ("outDir,O", value<path>(&out_dir) -> default_value(current_path().string()), "output directory")
            ("outBasename,o", value<string>(&out_basename), "specify the basename for output file(s)")
            ("cleanFastq,F", value< vector<path> >(&clean_fq) -> multitoken(), "cleaned fastq file name(s)")
        ;

        variables_map vm; 
        store(command_line_parser(argc, argv).options(desc).run(), vm);
        inter_dependency(vm, "outBasename", "outDir");
        inter_independency(vm, "cleanFastq", "outBasename");
        notify(vm);    

        if (vm.count("help")) {
            cout << desc << "\n";
            return 0;
        }

        for (int i = 0; i < raw_fq.size(); i++) {
            try {
                raw_fq[i] = canonical(raw_fq[i]);
            }
            catch (filesystem_error& e) {
                cerr << "error: No such file or directory: " << e.path1() << endl;
                return 1;
            }
        }
        
        try {
            out_dir = canonical(out_dir);
        }
        catch (filesystem_error& e) {
            cerr << "error: No such file or directory: " << e.path1() << endl;
            return 1;
        }

        if (vm.count("out-basename")) {
            for (int i = 0; i < raw_fq.size(); i++) {
                clean_fq.push_back(out_dir / path(out_basename + "_" + to_string(i + 1) + ".fastq.gz"));
            }
        }
        else {
            for (int i = 0; i < clean_fq.size(); i++) {
                clean_fq[i] = canonical(clean_fq[i].parent_path()) / clean_fq[i].filename();
            }
        }

        if (raw_fq.size() != clean_fq.size()) {
            cerr << "error: unequal numbers of raw and clean fastq: "
                << raw_fq.size() << " raw fastq and "
                << clean_fq.size() << " clean fastq." << endl;
            return 1;
        }
        
        cout << log_title() << "Welcome to filterfq!" << endl
            << log_title() << "INFO Command arguments -- filterfq ";
        for (variables_map::iterator i = vm.begin(); i != vm.end(); i++) {
            cout << i -> first << endl;
        }
        cout << log_title() << "INFO Filtering parameters -- Maximum rate of N bases: ";
        cout << "N base rate: " << base_N_rate << endl;
        cout << "average quality: " << ave_quality << endl;
        cout << "per base quality: " << base_quality << endl;
        cout << "low quality rate: " << low_quality_rate << endl;
        cout << "raw fastq: ";
        for (path p : raw_fq) {
            cout << p << " ";
        }
        cout << endl;
        cout << endl << "Output: " << endl;
        cout << "out dir: " << out_dir << endl;
        cout << "cleaned fastq quality system: " << clean_quality_sys << endl; 
        cout << "clean fastq: ";
        for (path p : clean_fq) {
            cout << p << " ";
        }
        cout << endl;

        // cout << "debug: " << endl;

        if (!prefer_specified_raw_quality_sys) {
            raw_quality_sys = check_quality_system(raw_fq[0]);
        }
        cout << "quality system: " << quality_sys[raw_quality_sys] << endl;

        // switch(raw_fq.size()) {
        //     case 1: 
        //         {
        //             ifstream infq(raw_fq[0].c_str(), ios_base::in | ios_base::binary);
        //             filtering_istream infq_decompressor;
        //             infq_decompressor.push(gzip_decompressor());
        //             infq_decompressor.push(infq);

        //             ofstream outfq(clean_fq[0].c_str(), ios_base::out | ios_base::binary);
        //             filtering_ostream outfq_compressor;
        //             outfq_compressor.push(gzip_compressor());
        //             outfq_compressor.push(outfq);

        //             string read_id_line, read_line, plus_line, quality_line;
        //             while (getline(infq_decompressor, read_id_line)) {
        //                 getline(infq_decompressor, read_line);
        //                 getline(infq_decompressor, plus_line);
        //                 getline(infq_decompressor, quality_line);
        //                 
        //                 cout << read_id_line << endl
        //                     << read_line << endl 
        //                     << plus_line << endl
        //                     << quality_line << endl;
        //                 cout << get_base_N_rate(read_line) << " "
        //                     << get_average_quality(quality_line, raw_quality_sys) << " "
        //                     << get_low_quality_rate(quality_line, raw_quality_sys, base_quality) << endl;
        //                 if (get_base_N_rate(read_line) > base_N_rate || get_average_quality(quality_line, raw_quality_sys) < ave_quality || get_low_quality_rate(quality_line, raw_quality_sys, base_quality) > low_quality_rate) {
        //                     cout << "delete!!" << endl;
        //                     continue;
        //                 }
        //                 else {
        //                     cout << "pass!!" << endl;
        //                     quality_system_convert(quality_line, raw_quality_sys, clean_quality_sys);
        //                     // cout << get_base_N_rate(read_line) << " "
        //                     //     << get_average_quality(quality_line, clean_quality_sys) << " "
        //                     //     << get_low_quality_rate(quality_line, clean_quality_sys, base_quality) << endl;
        //                     outfq_compressor << read_id_line << endl
        //                         << read_line << endl
        //                         << plus_line << endl
        //                         << quality_line << endl;
        //                 }
        //             }
        //             break;
        //         }
        //     case 2:
        //         {
        //             // iostream for the first fastq
        //             ifstream infq1(raw_fq[0].c_str(), ios_base::in | ios_base::binary);
        //             filtering_istream infq1_decompressor;
        //             infq1_decompressor.push(gzip_decompressor());
        //             infq1_decompressor.push(infq1);

        //             ofstream outfq1(clean_fq[0].c_str(), ios_base::out | ios_base::binary);
        //             filtering_ostream outfq1_compressor;
        //             outfq1_compressor.push(gzip_compressor());
        //             outfq1_compressor.push(outfq1);

        //             // iostream for the second fastq
        //             ifstream infq2(raw_fq[1].c_str(), ios_base::in | ios_base::binary);
        //             filtering_istream infq2_decompressor;
        //             infq2_decompressor.push(gzip_decompressor());
        //             infq2_decompressor.push(infq2);

        //             ofstream outfq2(clean_fq[1].c_str(), ios_base::out | ios_base::binary);
        //             filtering_ostream outfq2_compressor;
        //             outfq2_compressor.push(gzip_compressor());
        //             outfq2_compressor.push(outfq2);

        //             string read_id_line1, read_line1, plus_line1, quality_line1,
        //                    read_id_line2, read_line2, plus_line2, quality_line2;
        //             while (getline(infq1_decompressor, read_id_line1) && getline(infq2_decompressor, read_id_line2)) {
        //                 getline(infq1_decompressor, read_line1);
        //                 getline(infq1_decompressor, plus_line1);
        //                 getline(infq1_decompressor, quality_line1);
        //                 getline(infq2_decompressor, read_line2);
        //                 getline(infq2_decompressor, plus_line2);
        //                 getline(infq2_decompressor, quality_line2);

        //                 if (get_base_N_rate(read_line1) > base_N_rate || get_base_N_rate(read_line2) > base_N_rate 
        //                         || get_average_quality(quality_line1, raw_quality_sys) < ave_quality || get_average_quality(quality_line2, raw_quality_sys) < ave_quality 
        //                         || get_low_quality_rate(quality_line1, raw_quality_sys, base_quality) > low_quality_rate || get_low_quality_rate(quality_line2, raw_quality_sys, base_quality) > low_quality_rate) {
        //                     continue;
        //                 }
        //                 else {
        //                     quality_system_convert(quality_line1, raw_quality_sys, clean_quality_sys);
        //                     quality_system_convert(quality_line2, raw_quality_sys, clean_quality_sys);
        //                     outfq1_compressor << read_id_line1 << endl
        //                         << read_line1 << endl
        //                         << plus_line1 << endl
        //                         << quality_line1 << endl;
        //                     outfq2_compressor << read_id_line2 << endl
        //                         << read_line2 << endl
        //                         << plus_line2 << endl
        //                         << quality_line2 << endl;
        //                 }
        //             }
        //             break;
        //         }
        // }
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
