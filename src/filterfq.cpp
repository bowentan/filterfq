#include <boost/algorithm/string/join.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/thread.hpp>
#include <iostream>
#include <fstream>
#include <iterator>
#include <unordered_set>
#include <sstream>

#include <command_options.hpp>
#include <fastq_filter.hpp>
#include <quality_system.hpp>
#include <version.hpp>

// #define TESTING

using namespace boost::algorithm;
using namespace boost::program_options;
using namespace boost::iostreams;
using namespace boost::filesystem;
using namespace boost::posix_time;
using namespace boost::gregorian;
using namespace command_options;
using namespace fastq_filter;
using namespace quality_system;
using namespace std;

int main(int argc, char* argv[]) {
    try {
        ptime start_time = second_clock::local_time();
        // input variables
        const string quality_sys[5] = {"Sanger", "Solexa", "Illumina 1.3+", "Illumina 1.5+", "Illumina 1.8+"};
        const int THREAD_BLOCK_SIZE = 500000;
        bool only_get_read_info;
        bool prefer_specified_raw_quality_sys;
        // bool verbose;
        int n_thread;
        int max_read_len;
        int raw_quality_sys;
        int min_base_quality;
        float max_base_N_rate;
        float min_ave_quality;
        float max_low_quality_rate;
        path tmp_dir;
        vector<path> adapter;
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
            ("version,v", "print current version")
            ("thread,t", value<int>(&n_thread) -> default_value(8), "specify the number of threads to use, maximum 8")
            ("tmpDir,T", value<path>(&tmp_dir), "specify the directory to store temporary files, default as the same as \'outDir\'")
        ;
        
        options_description param("Input parameters & files", options_description::m_default_line_length * 2, options_description::m_default_line_length);
        param.add_options()
            ("checkQualitySystem,c", bool_switch(&only_get_read_info), "only check quality system of the fastq file")
            // ("verbose,v", bool_switch(&verbose), "print filtering information")
            ("maxReadLen,l", value<int>(&max_read_len) -> default_value(100), "maximum read length in the fastq file")
            ("baseNrate,N", value<float>(&max_base_N_rate) -> default_value(0.05), "maximum rate of \'N\' base allowed along a read")
            ("averageQuality,Q", value<float>(&min_ave_quality) -> default_value(0), "minimum average quality allowed along a read")
            ("perBaseQuality,q", value<int>(&min_base_quality) -> default_value(5), "minimum quality per base allowed along a read")
            ("lowQualityRate,r", value<float>(&max_low_quality_rate) -> default_value(0.5), "maximum low quality rate along a read")
            ("rawQualitySystem,s", value<int>(&raw_quality_sys), "specify quality system of raw fastq\n  0: Sanger\n  1: Solexa\n  2: Illumina 1.3+\n  3: Illumina 1.5+\n  4: Illumina 1.8+")
            ("preferSpecifiedRawQualitySystem,p", bool_switch(&prefer_specified_raw_quality_sys), "indicate that user prefers the given quality system to process")
            ("adapter,a", value< vector<path> >(&adapter) -> multitoken(), "adapter file(s)")
            ("rawFastq,f", value< vector<path> >(&raw_fq) -> required() -> multitoken(), "raw fastq file(s) that need cleaned, required")
        ;

        options_description output("Output parameters & files", options_description::m_default_line_length * 2, options_description::m_default_line_length);
        output.add_options()
            ("cleanQualitySystem,S", value<int>(&clean_quality_sys) -> default_value(4), "specify quality system of cleaned fastq, the same as rawQualitySystem")
            ("outDir,O", value<path>(&out_dir), "specify output directory")
            ("outBasename,o", value<string>(&out_basename), "specify the basename for output file(s)")
            // ("cleanFastq,F", value< vector<path> >(&clean_fq) -> multitoken(), "cleaned fastq file name(s), not used if outDir or outBasename is specified")
            // ("droppedFastq,D", value< vector<path> >(&dropped_fq) -> multitoken(), "fastq file(s) containing reads that are filtered out")
        ;
        
        stringstream usage;
        usage << "usage: VERSION " << VERSION << endl
            << "filterfq [-c] -f <fastq_1> [<fastq_2>] [OPTIONS]";
        options_description desc(usage.str());
        desc.add(generic).add(param).add(output);

        variables_map vm; 
        store(command_line_parser(argc, argv).options(desc).run(), vm);
        if (vm.count("help") || argc == 1) {
            cout << setprecision(2) << desc << endl;
            cout << endl << "Enquiry or bug report: Bowen Tan via Email (notebowentan@gmail.com) or GitHub (https://github.com/bowentan/filterfq)." << endl;
            return 0;
        }
        else if (vm.count("version")) {
            stringstream version;
            version << "Current version: " << VERSION;
            cout << version.str() << endl;
            return 0;
        }
        check_option_dependency(2, vm, "rawFastq", "outBasename", "outDir");
        check_option_dependency(1, vm, "outBasename", "outDir");
        check_option_dependency(1, vm, "outDir", "outBasename");
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

        if (only_get_read_info) {
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
            int* read_info = get_read_info(raw_fq[0]);
            cout << log_title() << "INFO -- After checking " << ((read_info[3] < 4000000) ? read_info[3] : 4000000) << " reads, min quality code is \'" 
                << (char)read_info[0] 
                << "\' and max quality code is \'" 
                << (char)read_info[1] 
                << "\', the quality system is probably " 
                << quality_sys[read_info[2]] 
                << ". " 
                << "The maximum length of scanned reads is " 
                << read_info[4] 
                << "." << endl;
            ptime end_time = second_clock::local_time();
            time_duration dt = end_time - start_time;
            cout << log_title() << "INFO -- Process finished successfully! "
                << dt.total_seconds() << " seconds elapsed. Thank you for using filterfq!" << endl;
            delete read_info;
            return 0;
        }
        
        try {
            out_dir = canonical(out_dir);
        }
        catch (filesystem_error& e) {
            cerr << "error: No such file or directory: " << e.path1().string() << endl;
            return 1;
        }

        // if (vm.count("outBasename")) {
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
        // }
        // else {
        //     for (int i = 0; i < clean_fq.size(); i++) {
        //         try {
        //             clean_fq[i] = canonical(clean_fq[i].parent_path()) / clean_fq[i].filename();
        //             dropped_fq[i] = canonical(dropped_fq[i].parent_path()) / dropped_fq[i].filename();
        //         }
        //         catch (filesystem_error& e) {
        //             cerr << "error: No such file or directory: " << e.path1().string() << endl;
        //             return 1;
        //         }
        //     }
        // }

        if (!vm.count("tmpDir")) {
            // tmp_dir = clean_fq[0].parent_path();
            tmp_dir = out_dir;
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
        
        if (vm.count("adapter")) {
            if (adapter.size() != raw_fq.size()) {
                cerr << "error: unequal numbers of raw fastq and adapter list: "
                    << raw_fq.size() << " raw fastq and "
                    << adapter.size() << " adapter list." << endl;
                return 1;
            }
            else {
                for (vector<path>::iterator p = adapter.begin(); p != adapter.end(); p++) {
                    try {
                        *p = canonical(*p);
                    }
                    catch (filesystem_error& e) {
                        cerr << "error: No such file or directory: " << e.path1().string() << endl;
                    }
                }
            }
        }

        cout << log_title() << "Welcome to filterfq!" << endl
            << log_title() << "INFO -- Arguments: ";
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
            else if (v.value().type() == typeid(string)) {
                cout << i -> first << "=" << v.as<string>() << " ";
            }
            else if (v.value().type() == typeid(path)) {
                cout << i -> first << "=" << canonical(v.as<path>()).string() << " ";
            }
            else if (v.value().type() == typeid(vector<path>)) {
                vector<path> vv = v.as< vector<path> >();
                vector<string> ss;
                for (vector<path>::const_iterator p = vv.begin(); p != vv.end(); p++){
                    ss.push_back(canonical(*p).string());
                }
                cout << i -> first << "=" << join(ss, ",") << " ";
            }
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

        int* read_info = get_read_info(raw_fq[0]);
        int read_len = read_info[3];
        cout << log_title() << "INFO -- After checking " << ((read_info[3] < 4000000) ? read_info[3] : 4000000) << " reads, min quality code is \'" 
            << (char)read_info[0] 
            << "\' and max quality code is \'" 
            << (char)read_info[1] 
            << "\', the quality system is probably " 
            << quality_sys[read_info[2]] 
            << ". " 
            << "The maximum length of scanned reads is " 
            << read_info[4] 
            << "." << endl;
        if (prefer_specified_raw_quality_sys) {
            cout << log_title() << "WARN -- User prefered specified quality system "
                << quality_sys[raw_quality_sys]
                << ". The program will treat quality codes correspondingly." << endl;
        }
        else {
            raw_quality_sys = read_info[2];
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

        if (read_info[4] > max_read_len) {
            max_read_len = read_info[4];
            cout << log_title() << "WARN -- The maximum read length exceeds the given maximum read length (100), change it to " << max_read_len << "." << endl;
        }

        if (n_thread > 8) {
            cout << log_title() << "WARN -- The given number of threads exceeds the maximum (8), changed it to 8." << endl;
            n_thread = 8;
        }

        if (read_info[3] < THREAD_BLOCK_SIZE * n_thread) {
            cout << log_title() << "WARN -- " << n_thread << " threads are redundant for filtering the given fastq(s), it is automatically adjusted to ";
            n_thread = (read_info[3] / THREAD_BLOCK_SIZE == 0) ? read_info[3] / THREAD_BLOCK_SIZE + 1 : read_info[3] / THREAD_BLOCK_SIZE;
            cout << n_thread << " threads in accordance with the given fastq(s)." << endl;
        }
        delete read_info;

        cout << log_title() << "INFO -- Start filtering..." << endl;
#ifdef TESTING
        return 0;
#endif

        statistic counter(raw_fq.size(), max_read_len);
        int param_int[4] = {min_base_quality, raw_quality_sys, clean_quality_sys, max_read_len};
        float param_float[3] = {max_base_N_rate, min_ave_quality, max_low_quality_rate};
        boost::thread t[n_thread];
        int thread_info[n_thread][3];

        vector< unordered_set<string> > adapter_read_id_lists;
        if (adapter.size() != 0) {
            for (vector<path>::iterator p = adapter.begin(); p != adapter.end(); p++)
                adapter_read_id_lists.push_back(load_adapter(*p));
        }

        for (int i = 0; i < n_thread; i++) {
            thread_info[i][0] = THREAD_BLOCK_SIZE;
            thread_info[i][1] = n_thread;
            thread_info[i][2] = i;
            t[i] = boost::thread(processor, 
                    raw_fq, 
                    clean_fq, 
                    dropped_fq, 
                    tmp_dir,
                    adapter_read_id_lists, 
                    param_int, 
                    param_float, 
                    &counter, 
                    thread_info[i]);
        }

        for (int i = 0; i < n_thread; i++)
            t[i].join();

        write_statistic(counter, out_dir);

        cout << log_title() << "INFO -- Merging tmp files..." << endl;
        merge(clean_fq, dropped_fq, tmp_dir, n_thread);
        cout << log_title() << "INFO -- Merge completed!" << endl;

        ptime end_time = second_clock::local_time();
        time_duration dt = end_time - start_time;
        cout << log_title() << "INFO -- Process finished successfully! "
            << dt.total_seconds() << " seconds elapsed. Thank you for using filterfq!" << endl;
    }
    catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        return 1;
    } catch(...) {
        cerr << "Exception of unknown type!\n";
    }

    return 0;
}
