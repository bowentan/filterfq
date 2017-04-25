#include <boost/program_options.hpp>
#include <cstdarg>
#include <iostream>
#include <iomanip>

#include <version.hpp>

namespace command_options {
    void print_usage() {
        std::cout << std::endl;
        std::cout << "Name: filterfq (filter/clean fastq sequencing reads)" << std::endl;
        std::cout << "Version: " << VERSION << std::endl;
        std::cout << "Enquiry contact: Bowen Tan (notebowentan@gmail.com)" << std::endl;
        std::cout << "GitHub page: https://github.com/bowentan/filterfq" << std::endl;
        std::cout << std::endl;
        std::cout << "usage: filterfq [-c] -f <fastq_1> [<fastq_2>] [OPTIONS]" << std::endl;
        std::cout << std::endl;
        std::cout << "General options:" << std::endl;
        std::cout << std::setw(30) << std::left << "  -h, --help" << std::setw(12) << " " << std::left << "print help message" << std::endl;
        std::cout << std::setw(30) << std::left << "  -v, --version" << std::setw(12) << " " << std::left << "print current version" << std::endl;
        std::cout << std::setw(30) << std::left << "  -T, --tmpDir" << std::setw(12) << "[<outDir>]" << std::left << "specify the directory to store temporary files" << std::endl;
        std::cout << std::setw(30) << std::left << "  -t, --thread" << std::setw(12) << "[8]" << std::left << "specify the number of threads to use, maximum 8" << std::endl;
        std::cout << std::endl;
        std::cout << "Input options:" << std::endl;
        std::cout << std::setw(30) << std::left << "  -f, --rawFastq" << std::setw(12) << " " << std::left << "raw fastq file(s) that cleaned. Required" << std::endl;
        std::cout << std::setw(30) << std::left << "  -a, --adapter" << std::setw(12) << " " << std::left << "adapter file(s) corresponding to given fastq file(s)" << std::endl;
        std::cout << std::setw(30) << std::left << "  -c, --checkQualitySystem" << std::setw(12) << " " << std::left << "only check quality system of give fastq(s). When not" << std::endl;
        std::cout << std::setw(30) << std::left << "  " << std::setw(12) << " " << std::left << "specified, filterfq will automatically check quality" << std::endl;
        std::cout << std::setw(30) << std::left << "  " << std::setw(12) << " " << std::left << "system before filtering" << std::endl;
        std::cout << std::setw(30) << std::left << "  -s, --rawQualitySystem" << std::setw(12) << " " << std::left << "specify quality system of given fastq(s)" << std::endl;
        std::cout << std::setw(30) << " " << std::setw(12) << " " << "    0: Sanger"<< std::endl;
        std::cout << std::setw(30) << " " << std::setw(12) << " " << "    1: Solexa"<< std::endl;
        std::cout << std::setw(30) << " " << std::setw(12) << " " << "    2: Illumina 1.3+"<< std::endl;
        std::cout << std::setw(30) << " " << std::setw(12) << " " << "    4: Illumina 1.8+"<< std::endl;
        std::cout << std::setw(30) << std::left << "  -p, --preferRawQuality" << std::setw(12) << " " << std::left << "output filtered fastq with the same quality system" << std::endl;
        std::cout << std::setw(30) << std::left << " " << std::setw(12) << " " << std::left << "with input fastq(s)" << std::endl;
        std::cout << std::setw(30) << std::left << "  -N, --baseNrate" << std::setw(12) << "[0.05]" << std::left << "maximum rate of 'N' base along a read" << std::endl;
        std::cout << std::setw(30) << std::left << "  -Q, --averageQuality" << std::setw(12) << "[0]" << std::left << "minimum average quality along a read" << std::endl;
        std::cout << std::setw(30) << std::left << "  -q, --baseQuality" << std::setw(12) << "[5]" << std::left << "minimum quality per base along a read" << std::endl;
        std::cout << std::setw(30) << std::left << "  -r, --lowQualityRate" << std::setw(12) << "[0.5]" << std::left << "maximum low quality rate along a read" << std::endl;
        std::cout << std::setw(30) << std::left << "  -m, --trim" << std::setw(12) << "[0]" << std::left << "the number of bases that should be trimmed at both" << std::endl;
        std::cout << std::setw(30) << std::left << " " << std::setw(12) << " " << std::left << "ends of a read" << std::endl;
        std::cout << std::setw(30) << " " << std::setw(12) << " " << "    " << std::setw(11) << std::left << "m" << ": all ends of reads are trimmed with"<< std::endl;
        std::cout << std::setw(30) << " " << std::setw(12) << " " << "    " << std::setw(11) << std::left << " " << "  m bases"<< std::endl;
        std::cout << std::setw(30) << " " << std::setw(12) << " " << "    " << std::setw(11) << std::left << "m n" << ": for pairend reads, both ends of"<< std::endl;
        std::cout << std::setw(30) << " " << std::setw(12) << " " << "    " << std::setw(11) << std::left << " " << "  fastq_1 and fastq_2 reads are"<< std::endl;
        std::cout << std::setw(30) << " " << std::setw(12) << " " << "    " << std::setw(11) << std::left << " " << "  trimmed with m and n bases"<< std::endl;
        std::cout << std::setw(30) << " " << std::setw(12) << " " << "    " << std::setw(11) << std::left << "m+ n-" << ": for single end read, trim m and n" << std::endl;
        std::cout << std::setw(30) << " " << std::setw(12) << " " << "    " << std::setw(11) << std::left << " " << "  bases at 5'(+) and 3'(-) ends of" << std::endl;
        std::cout << std::setw(30) << " " << std::setw(12) << " " << "    " << std::setw(11) << std::left << " " << "  of read; for pair end reads, trim" << std::endl;
        std::cout << std::setw(30) << " " << std::setw(12) << " " << "    " << std::setw(11) << std::left << " " << "  m and n bases at 5'(+) and 3'(-)" << std::endl;
        std::cout << std::setw(30) << " " << std::setw(12) << " " << "    " << std::setw(11) << std::left << " " << "  ends of both reads" << std::endl;
        std::cout << std::setw(30) << " " << std::setw(12) << " " << "    " << std::setw(11) << std::left << "m p+ q-" << ": only for pair end reads, the value" << std::endl;
        std::cout << std::setw(30) << " " << std::setw(12) << " " << "    " << std::setw(11) << std::left << " " << "  with no sign is the number of bases" << std::endl;
        std::cout << std::setw(30) << " " << std::setw(12) << " " << "    " << std::setw(11) << std::left << " " << "  trimmed at both ends; please put" << std::endl;
        std::cout << std::setw(30) << " " << std::setw(12) << " " << "    " << std::setw(11) << std::left << " " << "  values with signs next to each other" << std::endl;
        std::cout << std::setw(30) << " " << std::setw(12) << " " << "    " << std::setw(11) << std::left << "m+ n- p+ q-" << ": specify trimmed bases of each end,"<< std::endl;
        std::cout << std::setw(30) << " " << std::setw(12) << " " << "    " << std::setw(11) << std::left << " " << "  first '+' is for 5' end of fastq_1"<< std::endl;
        std::cout << std::setw(30) << " " << std::setw(12) << " " << "    " << std::setw(11) << std::left << " " << "  read and first '-' is for 3' end of"<< std::endl;
        std::cout << std::setw(30) << " " << std::setw(12) << " " << "    " << std::setw(11) << std::left << " " << "  fastq_1 read"<< std::endl;
        std::cout << std::setw(30) << std::left << "  -l, --minReadLen" << std::setw(12) << "[90]" << std::left << "minimum read length in ouput fastq(s), for trimming" << std::endl;
        std::cout << std::setw(30) << std::left << "  -L, --maxReadLen" << std::setw(12) << "[100]" << std::left << "maximum read length in input fastq(s), for statistical use" << std::endl;
        std::cout << std::endl;
        std::cout << "Output options:" << std::endl;
        std::cout << std::setw(30) << std::left << "  -S, --cleanQualitySystem" << std::setw(12) << "[4]" << std::left << "specify quality system of cleaned fastq(s)" << std::endl;
        std::cout << std::setw(30) << std::left << "  -o, --outBasename" << std::setw(12) << " " << std::left << "basename for output files. Required when filtering" << std::endl;
        std::cout << std::setw(30) << std::left << "  -O, --outDir" << std::setw(12) << " " << std::left << "output directory. Required when filtering" << std::endl;
        std::cout << std::endl;
    }

    // up to 2
    void check_option_dependency(int n, const boost::program_options::variables_map& vm, const char* for_what, ...) {
        if (!vm.count("help") && !vm["checkQualitySystem"].as<bool>()) {
            if (vm.count(for_what) && !vm[for_what].defaulted()) {
                std::vector<const char*> required_options;
                va_list vl;
                va_start(vl, for_what);
                for (int i = 0; i < n; i++) {
                    required_options.push_back(va_arg(vl, const char*));
                }
                va_end(vl);

                bool depended_exist = false;
                for (std::vector<const char*>::const_iterator i = required_options.begin(); i != required_options.end(); i++) {
                    if (vm[*i].value().type() == typeid(bool)) {
                        if (vm[*i].as<bool>()) {
                            depended_exist = true;
                            break;
                        }
                    }
                    else {
                        if (vm[*i].defaulted() || vm.count(*i)) {
                            depended_exist = true;
                            break;
                        }
                    }
                }

                if (!depended_exist) {
                    std::string s = std::string("Option '") + for_what + "' requires the combination with ";
                    switch (n) {
                        case 1:
                            {
                                s += "'" + std::string(required_options[0]) + "'.";
                                break;
                            }
                        case 2:
                            {
                                s += "'" + std::string(required_options[0]) + "' or '" + std::string(required_options[1]) + "'.";
                                break;
                            }
                        default:
                            {
                                for (std::vector<const char*>::const_iterator i = required_options.begin(); i != required_options.end(); i++) {
                                    if (i + 1 != required_options.end()) {
                                        s += "'" + std::string(*i) + "', ";
                                    }
                                    else {
                                        s += "or '" + std::string(*i) + "'.";
                                    }
                                }
                            }
                    }
                    throw std::logic_error(s);
                }
            }
        }
    }

    void check_option_independency(int n, const boost::program_options::variables_map& vm, ...) {
        if (!vm.count("help")) {
            std::vector<const char*> ops;
            va_list vl;
            va_start(vl, vm);
            for (int i = 0; i < n; i++) {
                ops.push_back(va_arg(vl, const char*));
            }
            va_end(vl);

            bool default_exist = false;
            for (std::vector<const char*>::const_iterator i = ops.begin(); i != ops.end(); i++) {
                if (vm[*i].defaulted()) {
                    default_exist = true;
                    break;
                }
            }

            if (!default_exist) {
                bool all_counted = true;
                for (std::vector<const char*>::const_iterator i = ops.begin(); i != ops.end(); i++) {
                    if (!vm.count(*i)) {
                        all_counted = false;
                        break;
                    }
                }
                if (all_counted) {
                    std::string s = std::string("Conflicting options ");
                    if (n == 2) {
                        s += "'" + std::string(ops[0]) + "' and '" + std::string(ops[1]) + "'.";
                    }
                    else {
                        for (std::vector<const char*>::const_iterator i = ops.begin(); i != ops.end(); i++) {
                            if (i + 1 != ops.end()) {
                                s += "'" + std::string(*i) + "', ";
                            }
                            else {
                                s += "and '" + std::string(*i) + "'.";
                            }
                        }
                    }
                    throw std::logic_error(s);
                }
            }
        }
    }
    
    int * parse_trim_param(std::vector<std::string> & s, int n_end) {
        int * trim_num;
        if (s.size() == 1) {
            int temp = std::stoi(s[0]);
            trim_num = new int[n_end * 2];
            for (int i = 0; i < n_end * 2; i++) {
                trim_num[i] = temp;
            }
        }
        // only two values given
        else if (s.size() == 2) {
            // their ends are the same or anyone of them is specified without 5' or 3'
            if ((s[0].back() == '+' || s[0].back() == '-') 
                    && (s[1].back() == '+' || s[1].back() == '-')) {
                if (s[0].back() == s[1].back()) {
                    throw std::logic_error("unrecognized arguments for '--trim/-m'");
                }
                else {
                    if (n_end == 1) {
                        if (s[0].back() == '+') {
                            trim_num = new int[n_end * 2]{std::stoi(s[0]), std::stoi(s[1])};
                        }
                        else {
                            trim_num = new int[n_end * 2]{std::stoi(s[1]), std::stoi(s[0])};
                        }
                    }
                    else {
                        int trim_5end, trim_3end;
                        if (s[0].back() == '+') {
                            trim_5end = std::stoi(s[0]);
                            trim_3end = std::stoi(s[1]);
                        }
                        else {
                            trim_5end = std::stoi(s[1]);
                            trim_3end = std::stoi(s[0]);
                        }
                        trim_num = new int[n_end * 2]{trim_5end, trim_3end, trim_5end, trim_3end};
                    }
                }
            }
            else if (s[0].back() >= '0' && s[0].back() <= '9' 
                    && s[1].back() >= '0' && s[1].back() <= '9') {
                if (n_end == 1) {
                    throw std::logic_error("unrecognized arguments for '--trim/-m'");
                }
                else {
                    int trim_pair1 = std::stoi(s[0]);
                    int trim_pair2 = std::stoi(s[1]);
                    trim_num = new int[n_end * 2]{trim_pair1, trim_pair1, trim_pair2, trim_pair2};
                }
            }
            else {
                throw std::logic_error("unrecognized arguments for '--trim/-m'");
            }
        }
        else if (s.size() == 3) {
            int nosign = 0;
            for (int i = 0; i < s.size(); i++) {
                if (s[i].back() != '+' && s[i].back() != '-') {
                    nosign++;
                }
            }
            if (nosign != 1 || (s[1].back() != '+' && s[1].back() != '-')) {
                throw std::logic_error("unrecognized arguments for '--trim/-m'");
            }

            if (s[0].back() != '+' && s[0].back() != '-') {
                int trim_pair1 = std::stoi(s[0]);
                int trim_pair2_5end, trim_pair2_3end;
                if (s[1].back() == '+') {
                    if (s[2].back() != '-') {
                        throw std::logic_error("unrecognized arguments for '--trim/-m'");
                    }
                    else {
                        trim_pair2_5end = std::stoi(s[1]);
                        trim_pair2_3end = std::stoi(s[2]);
                    }
                }
                else {
                    if (s[2].back() != '+') {
                        throw std::logic_error("unrecognized arguments for '--trim/-m'");
                    }
                    else {
                        trim_pair2_5end = std::stoi(s[2]);
                        trim_pair2_3end = std::stoi(s[1]);
                    }
                }
                trim_num = new int[n_end * 2]{trim_pair1, trim_pair1, trim_pair2_5end, trim_pair2_3end};
            }
            else if (s[2].back() != '+' && s[2].back() != '-') {
                int trim_pair2 = std::stoi(s[2]);
                int trim_pair1_5end, trim_pair1_3end;
                if (s[0].back() == '+') {
                    trim_pair1_5end = std::stoi(s[0]);
                    trim_pair1_3end = std::stoi(s[1]);
                }
                else {
                    trim_pair1_5end = std::stoi(s[1]);
                    trim_pair1_3end = std::stoi(s[0]);
                }
                trim_num = new int[n_end * 2]{trim_pair1_5end, trim_pair1_3end, trim_pair2, trim_pair2};
            }
        }
        else {
            int plus_sign = 0;
            int minus_sign = 0;
            for (int i = 0; i < n_end * 2; i++) {
                if (s[i].back() == '+') {
                    plus_sign++;
                }
                if (s[i].back() == '-') {
                    minus_sign++;
                }
            }
            if (plus_sign != 2 || minus_sign != 2) {
                throw std::logic_error("unrecognized arguments for '--trim/-m'");
            }
            int trim_5end[2] = {-1, -1};
            int trim_3end[2] = {-1, -1};
            for (int i = 0; i < n_end * 2; i++) {
                if (s[i].back() == '+') {
                    if (trim_5end[0] == -1) {
                        trim_5end[0] = std::stoi(s[i]);
                    }
                    else {
                        trim_5end[1] = std::stoi(s[i]);
                    }
                }
                else if (s[i].back() == '-') {
                    if (trim_3end[0] == -1) {
                        trim_3end[0] = std::stoi(s[i]);
                    }
                    else {
                        trim_3end[1] = std::stoi(s[i]);
                    }
                }
            }
            trim_num = new int[n_end * 2]{trim_5end[0], trim_3end[0], trim_5end[1], trim_3end[1]};
        }

        return trim_num;
    }
}
