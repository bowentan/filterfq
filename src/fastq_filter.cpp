#include <iostream>
#include <fstream>
#include <unordered_set>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/algorithm/minmax_element.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/thread.hpp>
#include <quality_system.hpp>

boost::mutex mutex;

namespace fastq_filter {

    struct statistic {
        unsigned long n_total;
        unsigned long n_filtered;
        unsigned long n_clean;
        std::vector< std::vector<unsigned long> > read_len_info;
        std::vector< std::vector<unsigned long> > filtered_read_info;
        std::vector< std::vector< std::vector<unsigned long> > > base_info;
        std::vector< std::vector< std::vector<unsigned long> > > base_quality_info;

        statistic(int n_end, int max_read_len) {
            n_filtered = 0;
            n_total = 0;
            n_clean = 0;
            read_len_info = std::vector< std::vector<unsigned long> >(2 * n_end, std::vector<unsigned long>(max_read_len));
            filtered_read_info = std::vector< std::vector<unsigned long> >(n_end, std::vector<unsigned long>(5));
            base_info = std::vector< std::vector< std::vector<unsigned long> > >(n_end, std::vector< std::vector<unsigned long> >(max_read_len, std::vector<unsigned long>(10)));
            base_quality_info = std::vector< std::vector< std::vector<unsigned long> > >(n_end * 2, std::vector< std::vector<unsigned long> >(max_read_len, std::vector<unsigned long>(42)));
        }
    };

    std::string log_title() {return "[filterfq | " + to_simple_string(boost::posix_time::second_clock::local_time()) + "] ";}

    int* check_quality_system(const boost::filesystem::path& filepath) {
        int* results = new int[5];
        std::ifstream infile (filepath.string(), std::ios_base::in | std::ios_base::binary);
        boost::iostreams::filtering_istream decompressor;
        decompressor.push(boost::iostreams::gzip_decompressor());
        decompressor.push(infile);
        std::string line;

        unsigned char min = '~';
        unsigned char max = '!';
        int n_read = 0;
        int max_len = 0;
        std::pair<std::string::iterator, std::string::iterator> extremum;
        // boost::random::mt19937 gen;
        // boost::random::uniform_int_distribution<> dist(1, 10);
        while (n_read <= 4000000 && getline(decompressor, line)) {
            for (int i = 0; i < 3; i++)
                getline(decompressor, line);
            // if (dist(gen) >= 5) 
            //     continue;

            if (line.length() > max_len) {
                max_len = line.length();
            }

            extremum = boost::minmax_element(line.begin(), line.end());
            if (*extremum.second >= max) {
                max = *extremum.second;
            }
            if (*extremum.first <= min) {
                min = *extremum.first;
            }
            n_read++;
        }
        results[0] = min;
        results[1] = max;
        boost::iostreams::close(decompressor, std::ios_base::in);
        if (min < ';') {
            if (max == 'I') {
                // sanger
                results[2] = 0;
            }
            else {
                // prefer illumina 1.8+
                results[2] = 4;
            }
        }
        else {
            if (min < '@') {
                // Solexa
                results[2] = 1;
            }
            else if (min < 'B') {
                // illumina 1.3+
                results[2] = 2;
            }
            else {
                // prefer illumina 1.5+
                results[2] = 3;
            }
        }
        results[3] = n_read;
        results[4] = max_len;
        return results;
    }
    
    float get_base_N_rate(const std::string& read_seq) {
        return std::count(read_seq.begin(), read_seq.end(), 'N') * 1.0 / read_seq.length(); }

    float get_average_quality(const std::string& quality_seq, const int sys) {
        unsigned int sum = 0;
        for (std::string::const_iterator s = quality_seq.begin(); s != quality_seq.end(); s++) {
            sum += *s;
        }
        return sum * 1.0 / quality_seq.size() - (quality_system::zero_quality)[sys];
        // corresponding quality min
        // boost::accumulators::accumulator_set<int, boost::accumulators::stats<boost::accumulators::tag::mean> > acc;
        // for(char q : quality_seq) {
        //     acc(q);
        // }
        // return boost::accumulators::mean(acc) - (quality_system::zero_quality)[sys];
    }

    float get_low_quality_rate(const std::string& quality_seq, const int sys, const int base_quality_threshold) {
        return std::count_if(quality_seq.begin(), quality_seq.end(), [sys, base_quality_threshold](char q) {return q - (quality_system::zero_quality)[sys] < base_quality_threshold;}) * 1.0 / quality_seq.length();
    }

    int* get_base_info(const std::string& base_seq) {
        int read_len = base_seq.length();
        int* base_info = new int[read_len + 1];
        base_info[0] = read_len;
        for (int i = 1; i < read_len + 1; i++) {
           switch(base_seq[i - 1]) {
               case 'A':
                   base_info[i] = 0;
                   break;
               case 'C':
                   base_info[i] = 1;
                   break;
               case 'G':
                   base_info[i] = 2;
                   break;
               case 'T':
                   base_info[i] = 3;
                   break;
               case 'N':
                   base_info[i] = 4;
                   break;
           }
        }
        return base_info;
    }

    int* get_base_quality_info(const std::string& quality_seq, const int quality_sys) {
        int read_len = quality_seq.length();
        int* base_quality_info = new int[read_len + 1];
        base_quality_info[0] = read_len;
        for (int i = 1; i < read_len + 1; i++) {
            (quality_seq[i - 1] - quality_system::zero_quality[quality_sys] >= 0) ? base_quality_info[i] = quality_seq[i - 1] - quality_system::zero_quality[quality_sys] : base_quality_info[i] = 0;
        }
        return base_quality_info;
    }

    std::unordered_set<std::string> load_adapter(boost::filesystem::path& adapter_file) {
        std::unordered_set<std::string> adapter_read_id_list;

        std::ifstream inadapter(adapter_file.string(), std::ios_base::in | std::ios_base::binary);
        boost::iostreams::filtering_istream inadapter_decompressor;
        inadapter_decompressor.push(boost::iostreams::gzip_decompressor());
        inadapter_decompressor.push(inadapter);

        std::cout << log_title() << "INFO -- Loading adapter list(s)..." << std::endl;
        std::string line;
        std::vector<std::string> line_splited;
        getline(inadapter_decompressor, line); // read out the header 
        while (getline(inadapter_decompressor, line)) {
            boost::split(line_splited, line, boost::is_any_of("\t"));
            adapter_read_id_list.insert(line_splited[0]);
            if (adapter_read_id_list.size() % 10000 == 0) {
                std::cout << log_title() << "INFO -- "
                    << adapter_read_id_list.size() << " read ID have been loaded." << std::endl;
            }
        }
        std::cout << log_title() << "INFO -- Totally "
            << adapter_read_id_list.size() << " have been loaded." << std::endl;
        return adapter_read_id_list;
    }

    void processor(std::vector<boost::filesystem::path>& infiles,
            std::vector<boost::filesystem::path>& clean_outfiles,
            std::vector<boost::filesystem::path>& dropped_outfiles,
            boost::filesystem::path& tmp_dir,
            std::vector< std::unordered_set<std::string> >& adapter_read_id_lists,
            int* param_int,
            float* param_float,
            statistic* stat,
            int* thread_info) {
        int stepsize = thread_info[0];
        int n_thread = thread_info[1];
        int thread = thread_info[2];
        
        int min_base_quality = param_int[0];
        int raw_quality_sys = param_int[1];
        int clean_quality_sys = param_int[2];
        int max_read_len = param_int[3];

        float max_base_N_rate = param_float[0];
        float min_ave_quality = param_float[1];
        float max_low_quality_rate = param_float[2];

        if (infiles.size() == 1) {
            // file stream
            std::ifstream infq(infiles[0].string(), std::ios_base::in | std::ios_base::binary);
            boost::iostreams::filtering_istream infq_decompressor;
            infq_decompressor.push(boost::iostreams::gzip_decompressor());
            infq_decompressor.push(infq);

            std::ofstream clean_outfq((tmp_dir / clean_outfiles[0].filename()).string() + "." + std::to_string(thread) + ".tmp", std::ios_base::out | std::ios_base::binary); 
            boost::iostreams::filtering_ostream clean_outfq_compressor;
            clean_outfq_compressor.push(boost::iostreams::gzip_compressor());
            clean_outfq_compressor.push(clean_outfq);

            std::ofstream dropped_outfq((tmp_dir / dropped_outfiles[0].filename()).string() + "." + std::to_string(thread) + ".tmp", std::ios_base::out | std::ios_base::binary);
            boost::iostreams::filtering_ostream dropped_outfq_compressor;
            dropped_outfq_compressor.push(boost::iostreams::gzip_compressor());
            dropped_outfq_compressor.push(dropped_outfq);
            
            statistic local_counter(1, max_read_len);
            unsigned int step_counter = 0;
            int* base_info;
            int* base_quality_info;
            int* clean_base_quality_info;
            std::string read_id_line, read_line, plus_line, quality_line;
            bool is_filtered;

            for (int i = 0; i < stepsize * thread; i++) {
                for (int j = 0; j < 4; j++) {
                    getline(infq_decompressor, read_id_line);
                }
            }

            while (getline(infq_decompressor, read_id_line)) {
                getline(infq_decompressor, read_line);
                getline(infq_decompressor, plus_line);
                getline(infq_decompressor, quality_line);

                if (read_line.length() > max_read_len) {
                    std::cout << log_title() << "ERROR -- There are reads whose length (" << read_line.length() << ") exceeds the maximum read length (" << max_read_len << ") so that segmentation fault may occur. Please set the argument of the parameter \'-maxReadLen/-l\' as one integer larger than or equal to " << read_line.length() << "." << std::endl;
                    exit(1);
                }

                base_info = get_base_info(read_line);
                base_quality_info = get_base_quality_info(quality_line, raw_quality_sys);
                is_filtered = false;

                local_counter.read_len_info[0][base_info[0] - 1]++;
                if (get_base_N_rate(read_line) > max_base_N_rate) {
                    local_counter.filtered_read_info[0][0]++;
                    if (!is_filtered) {
                        local_counter.n_filtered++;
                        local_counter.filtered_read_info[0][4]++;
                        is_filtered = true;
                    }
                }
                if (get_average_quality(quality_line, raw_quality_sys) < min_ave_quality) {
                    local_counter.filtered_read_info[0][1]++;
                    if (!is_filtered) {
                        local_counter.n_filtered++;
                        local_counter.filtered_read_info[0][4]++;
                        is_filtered = true;
                    }
                }
                if (get_low_quality_rate(quality_line, raw_quality_sys, min_base_quality) > max_low_quality_rate) {
                    local_counter.filtered_read_info[0][2]++;
                    if (!is_filtered) {
                        local_counter.n_filtered++;
                        local_counter.filtered_read_info[0][4]++;
                        is_filtered = true;
                    }
                }
                if (adapter_read_id_lists.size() != 0) {
                    std::string s = read_id_line.substr(1, read_id_line.size() - 1);
                    if (adapter_read_id_lists[0].count(s) > 0) {
                        local_counter.filtered_read_info[0][3]++;
                        if (!is_filtered) {
                            local_counter.n_filtered++;
                            local_counter.filtered_read_info[0][4]++;
                            is_filtered = true;
                        }
                    }
                }

                local_counter.n_total++;
                for (int i = 1; i < base_info[0] + 1; i++) {
                    local_counter.base_info[0][i - 1][base_info[i]]++;
                    local_counter.base_quality_info[0][i - 1][base_quality_info[i]]++;
                }
                // if (verbose && counter[0] % 50000 == 0) {
                //     std::cout << log_title() << "INFO " 
                //         << std::fixed
                //         << std::setw(12) << std::setprecision(6) << counter[0] << " | "
                //         << std::setw(12) << std::setprecision(6) << counter[1] << " | "
                //         << std::setw(12) << std::setprecision(6) << counter[1] * 100.0 / counter[0] << " | "
                //         << std::setw(12) << std::setprecision(6) << counter[2] << " | "
                //         << std::setw(12) << std::setprecision(6) << counter[3] << " | "
                //         << std::setw(12) << std::setprecision(6) << counter[4] << " | "
                //         << std::setw(12) << std::setprecision(6) << counter[5] << " | "
                //         << std::endl;
                // }
                // mutex.unlock();

                if (!is_filtered) {
                    local_counter.n_clean++;
                    local_counter.read_len_info[1][base_info[0] - 1]++;
                    if (raw_quality_sys != clean_quality_sys) {
                        quality_system::quality_system_convert(quality_line, raw_quality_sys, clean_quality_sys);
                    }
                    clean_base_quality_info = get_base_quality_info(quality_line, clean_quality_sys);
                    for (int i = 1; i < clean_base_quality_info[0] + 1; i++) {
                        local_counter.base_info[0][i - 1][base_info[i] + 5]++;
                        local_counter.base_quality_info[1][i - 1][clean_base_quality_info[i]]++;
                    }
                    clean_outfq_compressor << read_id_line << std::endl
                        << read_line << std::endl
                        << plus_line << std::endl
                        << quality_line << std::endl;
                    delete [] clean_base_quality_info;
                }
                else {
                    if (raw_quality_sys != clean_quality_sys) {
                        quality_system::quality_system_convert(quality_line, raw_quality_sys, clean_quality_sys);
                    }
                    dropped_outfq_compressor << read_id_line << std::endl
                        << read_line << std::endl
                        << plus_line << std::endl
                        << quality_line << std::endl;
                }

                step_counter++;
                if (step_counter % stepsize == 0) {
                    for (int i = 0; i < stepsize * (n_thread - 1); i++) {
                        for (int j = 0; j < 4; j++) {
                            getline(infq_decompressor, read_id_line);
                        }
                    }
                }
                delete [] base_info;
                delete [] base_quality_info;
            }
            close(infq_decompressor, std::ios_base::in);
            close(clean_outfq_compressor, std::ios_base::out);
            close(dropped_outfq_compressor, std::ios_base::out);

            mutex.lock();
            (*stat).n_filtered += local_counter.n_filtered;
            (*stat).n_total += local_counter.n_total;
            (*stat).n_clean += local_counter.n_clean;
            for (int i = 0; i < infiles.size(); i++) {
                for (int j = 0; j < (*stat).read_len_info[i].size(); j++) {
                    (*stat).read_len_info[i][j] += local_counter.read_len_info[i][j];
                }
                for (int j = 0; j < (*stat).filtered_read_info[i].size(); j++) {
                    (*stat).filtered_read_info[i][j] += local_counter.filtered_read_info[i][j];
                }
                for (int j = 0; j < (*stat).base_info[i].size(); j++) {
                    for (int k = 0; k < (*stat).base_info[i][j].size(); k++) {
                        (*stat).base_info[i][j][k] += local_counter.base_info[i][j][k];
                    }
                    for (int k = 0; k < (*stat).base_quality_info[i][j].size(); k++) {
                        (*stat).base_quality_info[2 * i][j][k] += local_counter.base_quality_info[2 * i][j][k];
                        (*stat).base_quality_info[2 * i + 1][j][k] += local_counter.base_quality_info[2 * i + 1][j][k];
                    }
                }
            }
            mutex.unlock();
        }
        else if (infiles.size() == 2) {
            std::ifstream infq1(infiles[0].string(), std::ios_base::in | std::ios_base::binary);
            std::ifstream infq2(infiles[1].string(), std::ios_base::in | std::ios_base::binary);
            boost::iostreams::filtering_istream infq1_decompressor;
            boost::iostreams::filtering_istream infq2_decompressor;
            infq1_decompressor.push(boost::iostreams::gzip_decompressor());
            infq2_decompressor.push(boost::iostreams::gzip_decompressor());
            infq1_decompressor.push(infq1);
            infq2_decompressor.push(infq2);

            std::ofstream clean_outfq1((tmp_dir / clean_outfiles[0].filename()).string() + "." + std::to_string(thread) + ".tmp", std::ios_base::out | std::ios_base::binary); 
            std::ofstream clean_outfq2((tmp_dir / clean_outfiles[1].filename()).string() + "." + std::to_string(thread) + ".tmp", std::ios_base::out | std::ios_base::binary); 
            boost::iostreams::filtering_ostream clean_outfq1_compressor;
            boost::iostreams::filtering_ostream clean_outfq2_compressor;
            clean_outfq1_compressor.push(boost::iostreams::gzip_compressor());
            clean_outfq2_compressor.push(boost::iostreams::gzip_compressor());
            clean_outfq1_compressor.push(clean_outfq1);
            clean_outfq2_compressor.push(clean_outfq2);

            std::ofstream dropped_outfq1((tmp_dir / dropped_outfiles[0].filename()).string() + "." + std::to_string(thread) + ".tmp", std::ios_base::out | std::ios_base::binary);
            std::ofstream dropped_outfq2((tmp_dir / dropped_outfiles[1].filename()).string() + "." + std::to_string(thread) + ".tmp", std::ios_base::out | std::ios_base::binary);
            boost::iostreams::filtering_ostream dropped_outfq1_compressor;
            boost::iostreams::filtering_ostream dropped_outfq2_compressor;
            dropped_outfq1_compressor.push(boost::iostreams::gzip_compressor());
            dropped_outfq2_compressor.push(boost::iostreams::gzip_compressor());
            dropped_outfq1_compressor.push(dropped_outfq1);
            dropped_outfq2_compressor.push(dropped_outfq2);
            
            statistic local_counter(2, max_read_len);
            unsigned int step_counter = 0;
            int* base_info1;
            int* base_info2;
            int* base_quality_info1;
            int* base_quality_info2;
            int* clean_base_quality_info1;
            int* clean_base_quality_info2;
            std::string read_id_line1, read_line1, plus_line1, quality_line1;
            std::string read_id_line2, read_line2, plus_line2, quality_line2;
            bool is_filtered1;
            bool is_filtered2;
            bool is_pair_filtered;

            for (int i = 0; i < stepsize * thread; i++) {
                for (int j = 0; j < 4; j++) {
                    getline(infq1_decompressor, read_id_line1);
                    getline(infq2_decompressor, read_id_line2);
                }
            }

            while (getline(infq1_decompressor, read_id_line1) && getline(infq2_decompressor, read_id_line2)) {
                getline(infq1_decompressor, read_line1);
                getline(infq2_decompressor, read_line2);
                getline(infq1_decompressor, plus_line1);
                getline(infq2_decompressor, plus_line2);
                getline(infq1_decompressor, quality_line1);
                getline(infq2_decompressor, quality_line2);

                if (read_line1.length() > max_read_len || read_line2.length() > max_read_len) {
                    std::cout << log_title() << "ERROR -- There are reads whose length (" << std::max(read_line1.length(), read_line2.length()) << ") exceeds the maximum read length (" << max_read_len << ") so that segmentation fault may occur. Please set the argument of the parameter \'-maxReadLen/-l\' as one integer larger than or equal to " << std::max(read_line1.length(), read_line2.length()) << "." << std::endl;
                    exit(1);
                }

                base_info1 = get_base_info(read_line1);
                base_info2 = get_base_info(read_line2);
                base_quality_info1 = get_base_quality_info(quality_line1, raw_quality_sys);
                base_quality_info2 = get_base_quality_info(quality_line2, raw_quality_sys);
                is_filtered1 = false;
                is_filtered2 = false;
                is_pair_filtered = false;

                local_counter.read_len_info[0][base_info1[0] - 1]++;
                local_counter.read_len_info[2][base_info2[0] - 1]++;
                if (get_base_N_rate(read_line1) > max_base_N_rate) {
                    local_counter.filtered_read_info[0][0]++;
                    if (!is_filtered1) {
                        local_counter.filtered_read_info[0][4]++;
                        is_filtered1 = true;
                    }
                    if (!is_pair_filtered) {
                        local_counter.n_filtered++;
                        is_pair_filtered = true;
                    }
                }
                if (get_base_N_rate(read_line2) > max_base_N_rate) {
                    local_counter.filtered_read_info[1][0]++;
                    if (!is_filtered2) {
                        local_counter.filtered_read_info[1][4]++;
                        is_filtered2 = true;
                    }
                    if (!is_pair_filtered) {
                        local_counter.n_filtered++;
                        is_pair_filtered = true;
                    }
                }

                if (get_average_quality(quality_line1, raw_quality_sys) < min_ave_quality) {
                    local_counter.filtered_read_info[0][1]++;
                    if (!is_filtered1) {
                        local_counter.filtered_read_info[0][4]++;
                        is_filtered1 = true;
                    }
                    if (!is_pair_filtered) {
                        local_counter.n_filtered++;
                        is_pair_filtered = true;
                    }
                }
                if (get_average_quality(quality_line2, raw_quality_sys) < min_ave_quality) {
                    local_counter.filtered_read_info[1][1]++;
                    if (!is_filtered2) {
                        local_counter.filtered_read_info[1][4]++;
                        is_filtered2 = true;
                    }
                    if (!is_pair_filtered) {
                        local_counter.n_filtered++;
                        is_pair_filtered = true;
                    }
                }

                if (get_low_quality_rate(quality_line1, raw_quality_sys, min_base_quality) > max_low_quality_rate) {
                    local_counter.filtered_read_info[0][2]++;
                    if (!is_filtered1) {
                        local_counter.filtered_read_info[0][4]++;
                        is_filtered1 = true;
                    }
                    if (!is_pair_filtered) {
                        local_counter.n_filtered++;
                        is_pair_filtered = true;
                    }
                }
                if (get_low_quality_rate(quality_line2, raw_quality_sys, min_base_quality) > max_low_quality_rate) {
                    local_counter.filtered_read_info[1][2]++;
                    if (!is_filtered2) {
                        local_counter.filtered_read_info[1][4]++;
                        is_filtered2 = true;
                    }
                    if (!is_pair_filtered) {
                        local_counter.n_filtered++;
                        is_pair_filtered = true;
                    }
                }
                if (adapter_read_id_lists.size() != 0) {
                    std::string s1 = read_id_line1.substr(1, read_id_line1.size() - 1);
                    std::string s2 = read_id_line2.substr(1, read_id_line2.size() - 1);
                    if (adapter_read_id_lists[0].count(s1) > 0) {
                        local_counter.filtered_read_info[0][3]++;
                        if (!is_filtered1) {
                            local_counter.filtered_read_info[0][4]++;
                            is_filtered1 = true;
                        }
                        if (!is_pair_filtered) {
                            local_counter.n_filtered++;
                            is_pair_filtered = true;
                        }
                    }
                    if (adapter_read_id_lists[1].count(s2) > 0) {
                        local_counter.filtered_read_info[1][3]++;
                        if (!is_filtered2) {
                            local_counter.filtered_read_info[1][4]++;
                            is_filtered2 = true;
                        }
                        if (!is_pair_filtered) {
                            local_counter.n_filtered++;
                            is_pair_filtered = true;
                        }
                    }
                }

                local_counter.n_total++;
                for (int i = 1; i < base_info1[0] + 1; i++) {
                    local_counter.base_info[0][i - 1][base_info1[i]]++;
                    local_counter.base_quality_info[0][i - 1][base_quality_info1[i]]++;
                }
                for (int i = 1; i < base_info2[0] + 1; i++) {
                    local_counter.base_info[1][i - 1][base_info2[i]]++;
                    local_counter.base_quality_info[2][i - 1][base_quality_info2[i]]++;
                }
                // if (verbose && counter[0] % 50000 == 0) {
                //     std::cout << log_title() << "INFO " 
                //         << std::fixed
                //         << std::setw(12) << std::setprecision(6) << counter[0] << " | "
                //         << std::setw(12) << std::setprecision(6) << counter[1] << " | "
                //         << std::setw(12) << std::setprecision(6) << counter[1] * 100.0 / counter[0] << " | "
                //         << std::setw(12) << std::setprecision(6) << counter[2] << " | "
                //         << std::setw(12) << std::setprecision(6) << counter[3] << " | "
                //         << std::setw(12) << std::setprecision(6) << counter[4] << " | "
                //         << std::setw(12) << std::setprecision(6) << counter[5] << " | "
                //         << std::endl;
                // }
                // mutex.unlock();

                if (!is_pair_filtered) {
                    local_counter.n_clean++;
                    local_counter.read_len_info[1][base_info1[0] - 1]++;
                    local_counter.read_len_info[3][base_info2[0] - 1]++;
                    if (raw_quality_sys != clean_quality_sys) {
                        quality_system::quality_system_convert(quality_line1, raw_quality_sys, clean_quality_sys);
                        quality_system::quality_system_convert(quality_line2, raw_quality_sys, clean_quality_sys);
                    }
                    clean_base_quality_info1 = get_base_quality_info(quality_line1, clean_quality_sys);
                    clean_base_quality_info2 = get_base_quality_info(quality_line2, clean_quality_sys);
                    for (int i = 1; i < clean_base_quality_info1[0] + 1; i++) {
                        local_counter.base_info[0][i - 1][base_info1[i] + 5]++;
                        local_counter.base_quality_info[1][i - 1][clean_base_quality_info1[i]]++;
                    }
                    for (int i = 1; i < clean_base_quality_info2[0] + 1; i++) {
                        local_counter.base_info[1][i - 1][base_info2[i] + 5]++;
                        local_counter.base_quality_info[3][i - 1][clean_base_quality_info2[i]]++;
                    }
                    clean_outfq1_compressor << read_id_line1 << std::endl
                        << read_line1 << std::endl
                        << plus_line1 << std::endl
                        << quality_line1 << std::endl;
                    clean_outfq2_compressor << read_id_line2 << std::endl
                        << read_line2 << std::endl
                        << plus_line2 << std::endl
                        << quality_line2 << std::endl;
                    delete [] clean_base_quality_info1;
                    delete [] clean_base_quality_info2;
                }
                else {
                    if (raw_quality_sys != clean_quality_sys) {
                        quality_system::quality_system_convert(quality_line1, raw_quality_sys, clean_quality_sys);
                        quality_system::quality_system_convert(quality_line2, raw_quality_sys, clean_quality_sys);
                    }
                    dropped_outfq1_compressor << read_id_line1 << std::endl
                        << read_line1 << std::endl
                        << plus_line1 << std::endl
                        << quality_line1 << std::endl;
                    dropped_outfq2_compressor << read_id_line2 << std::endl
                        << read_line2 << std::endl
                        << plus_line2 << std::endl
                        << quality_line2 << std::endl;
                }
                delete [] base_info1;
                delete [] base_info2;
                delete [] base_quality_info1;
                delete [] base_quality_info2;

                step_counter++;
                if (step_counter % stepsize == 0) {
                    for (int i = 0; i < stepsize * (n_thread - 1); i++) {
                        for (int j = 0; j < 4; j++) {
                            getline(infq1_decompressor, read_id_line1);
                            getline(infq2_decompressor, read_id_line2);
                        }
                    }
                }
            }
            close(infq1_decompressor, std::ios_base::in);
            close(infq2_decompressor, std::ios_base::in);
            close(clean_outfq1_compressor, std::ios_base::out);
            close(clean_outfq2_compressor, std::ios_base::out);
            close(dropped_outfq1_compressor, std::ios_base::out);
            close(dropped_outfq2_compressor, std::ios_base::out);

            mutex.lock();
            (*stat).n_filtered += local_counter.n_filtered;
            (*stat).n_total += local_counter.n_total;
            (*stat).n_clean += local_counter.n_clean;
            for (int i = 0; i < infiles.size(); i++) {
                for (int j = 0; j < (*stat).read_len_info[i].size(); j++) {
                    (*stat).read_len_info[i][j] += local_counter.read_len_info[i][j];
                }
                for (int j = 0; j < (*stat).filtered_read_info[i].size(); j++) {
                    (*stat).filtered_read_info[i][j] += local_counter.filtered_read_info[i][j];
                }
                for (int j = 0; j < (*stat).base_info[i].size(); j++) {
                    for (int k = 0; k < (*stat).base_info[i][j].size(); k++) {
                        (*stat).base_info[i][j][k] += local_counter.base_info[i][j][k];
                    }
                    for (int k = 0; k < (*stat).base_quality_info[i][j].size(); k++) {
                        (*stat).base_quality_info[2 * i][j][k] += local_counter.base_quality_info[2 * i][j][k];
                        (*stat).base_quality_info[2 * i + 1][j][k] += local_counter.base_quality_info[2 * i + 1][j][k];
                    }
                }
            }
            mutex.unlock();
        }
    }

    void merge_single(boost::filesystem::path& outfile, boost::filesystem::path& tmp_dir, int n_thread) {
        boost::filesystem::remove(outfile);
        std::ofstream out(outfile.string(), std::ios_base::out | std::ios_base::binary | std::ios_base::app);
        
        for (int j = 0; j < n_thread; j++) {
            std::string tmp_filename = (tmp_dir / outfile.filename()).string() + "." + std::to_string(j) + ".tmp";
            std::ifstream tmp_file(tmp_filename, std::ios_base::in | std::ios_base::binary);

            out << tmp_file.rdbuf();

            tmp_file.close();

            boost::filesystem::remove(tmp_filename);
        }
        
        out.close();
    }

    void merge(std::vector<boost::filesystem::path>& clean_outfiles,
            std::vector<boost::filesystem::path>& dropped_outfiles,
            boost::filesystem::path& tmp_dir,
            int n_thread) {
        boost::thread t[clean_outfiles.size() * 2];
        for (int i = 0; i < clean_outfiles.size(); i++) {
            t[2 * i] = boost::thread(merge_single, clean_outfiles[i], tmp_dir, n_thread);
            t[2 * i + 1] = boost::thread(merge_single, dropped_outfiles[i], tmp_dir, n_thread);
        }

        for (int i = 0; i < clean_outfiles.size() * 2; i++) {
            t[i].join();
        }
    }

    void write_statistic(statistic stat, boost::filesystem::path& out_dir) {
        int n_end = stat.base_info.size();
        // filtered info
        std::ofstream filtered_read_info_file(out_dir.string() + "/Statistics_of_reads.txt", std::ios_base::out);
        if (n_end == 1) {
            filtered_read_info_file << "item\tfastq\t%_total\t%_actual_filtered\t%_marked_filtered" << std::endl;
            filtered_read_info_file << std::fixed << std::setprecision(2) << "total\t" << stat.n_total << '\t' << stat.n_total * 100.0 / stat.n_total << "\t-\t-" << std::endl;
            filtered_read_info_file << std::fixed << std::setprecision(2) << "actual_filtered\t" << stat.n_filtered << '\t' << stat.n_filtered * 100.0 / stat.n_total << '\t' << stat.n_filtered * 100.0 / stat.n_filtered << "\t-" << std::endl;
            filtered_read_info_file << std::fixed << std::setprecision(2) << "marked_filtered" 
                << '\t' << stat.filtered_read_info[0][4] 
                << '\t' << stat.filtered_read_info[0][4] * 100.0 / stat.n_total
                << '\t' << stat.filtered_read_info[0][4] * 100.0 / stat.n_filtered
                << '\t' << stat.filtered_read_info[0][4] * 100.0 / stat.filtered_read_info[0][4]
                << std::endl;
            filtered_read_info_file << std::fixed << std::setprecision(2) << "clean\t" << stat.n_clean << '\t' << stat.n_clean * 100.0 / stat.n_total << "\t-\t-" << std::endl;
            filtered_read_info_file << std::fixed << std::setprecision(2) << "high_N_rate"
                << '\t' << stat.filtered_read_info[0][0]
                << '\t' << stat.filtered_read_info[0][0] * 100.0 / stat.n_total
                << '\t' << stat.filtered_read_info[0][0] * 100.0 / stat.n_filtered
                << '\t' << stat.filtered_read_info[0][0] * 100.0 / stat.filtered_read_info[0][4]
                << std::endl;
            filtered_read_info_file << std::fixed << std::setprecision(2) << "low_ave_quality"
                << '\t' << stat.filtered_read_info[0][1]
                << '\t' << stat.filtered_read_info[0][1] * 100.0 / stat.n_total
                << '\t' << stat.filtered_read_info[0][1] * 100.0 / stat.n_filtered
                << '\t' << stat.filtered_read_info[0][1] * 100.0 / stat.filtered_read_info[0][4]
                << std::endl;
            filtered_read_info_file << std::fixed << std::setprecision(2) << "high_low_quality_rate"
                << '\t' << stat.filtered_read_info[0][2]
                << '\t' << stat.filtered_read_info[0][2] * 100.0 / stat.n_total
                << '\t' << stat.filtered_read_info[0][2] * 100.0 / stat.n_filtered
                << '\t' << stat.filtered_read_info[0][2] * 100.0 / stat.filtered_read_info[0][4]
                << std::endl;
            filtered_read_info_file << std::fixed << std::setprecision(2) << "with_adapter"
                << '\t' << stat.filtered_read_info[0][3]
                << '\t' << stat.filtered_read_info[0][3] * 100.0 / stat.n_total
                << '\t' << stat.filtered_read_info[0][3] * 100.0 / stat.n_filtered
                << '\t' << stat.filtered_read_info[0][3] * 100.0 / stat.filtered_read_info[0][4]
                << std::endl;
        }
        else if (n_end  == 2) {
            filtered_read_info_file << "item\tfastq_1\t%_total\t%_actual_filtered\t%_marked_filtered\tfastq_2\t%_total\t%_actual_filtered\t%_marked_filtered" << std::endl;
            filtered_read_info_file << std::fixed << std::setprecision(2) << "total\t" << stat.n_total << '\t' << stat.n_total * 100.0 / stat.n_total << "\t-\t-\t" << stat.n_total << '\t' << stat.n_total * 100.0 / stat.n_total << "\t-\t-" << std::endl;
            filtered_read_info_file << std::fixed << std::setprecision(2) << "actual_filtered\t" << stat.n_filtered << '\t' << stat.n_filtered * 100.0 / stat.n_total << '\t' << stat.n_filtered * 100.0 / stat.n_filtered << "\t-\t" << stat.n_filtered << '\t' << stat.n_filtered * 100.0 / stat.n_total << '\t' << stat.n_filtered * 100.0 / stat.n_filtered << "\t-" << std::endl;
            filtered_read_info_file << std::fixed << std::setprecision(2) << "marked_filtered" 
                << '\t' << stat.filtered_read_info[0][4] 
                << '\t' << stat.filtered_read_info[0][4] * 100.0 / stat.n_total
                << '\t' << stat.filtered_read_info[0][4] * 100.0 / stat.n_filtered
                << '\t' << stat.filtered_read_info[0][4] * 100.0 / stat.filtered_read_info[0][4]
                << '\t' << stat.filtered_read_info[1][4] 
                << '\t' << stat.filtered_read_info[1][4] * 100.0 / stat.n_total
                << '\t' << stat.filtered_read_info[1][4] * 100.0 / stat.n_filtered
                << '\t' << stat.filtered_read_info[1][4] * 100.0 / stat.filtered_read_info[1][4]
                << std::endl;
            filtered_read_info_file << std::fixed << std::setprecision(2) << "clean\t" << stat.n_clean << '\t' << stat.n_clean * 100.0 / stat.n_total << "\t-\t-\t" << stat.n_clean << '\t' << stat.n_clean * 100.0 / stat.n_total << "\t-\t-" << std::endl;
            filtered_read_info_file << std::fixed << std::setprecision(2) << "high_N_rate"
                << '\t' << stat.filtered_read_info[0][0]
                << '\t' << stat.filtered_read_info[0][0] * 100.0 / stat.n_total
                << '\t' << stat.filtered_read_info[0][0] * 100.0 / stat.n_filtered
                << '\t' << stat.filtered_read_info[0][0] * 100.0 / stat.filtered_read_info[0][4]
                << '\t' << stat.filtered_read_info[1][0] 
                << '\t' << stat.filtered_read_info[1][0] * 100.0 / stat.n_total
                << '\t' << stat.filtered_read_info[1][0] * 100.0 / stat.n_filtered
                << '\t' << stat.filtered_read_info[1][0] * 100.0 / stat.filtered_read_info[1][4]
                << std::endl;
            filtered_read_info_file << std::fixed << std::setprecision(2) << "low_ave_quality"
                << '\t' << stat.filtered_read_info[0][1]
                << '\t' << stat.filtered_read_info[0][1] * 100.0 / stat.n_total
                << '\t' << stat.filtered_read_info[0][1] * 100.0 / stat.n_filtered
                << '\t' << stat.filtered_read_info[0][1] * 100.0 / stat.filtered_read_info[0][4]
                << '\t' << stat.filtered_read_info[1][1] 
                << '\t' << stat.filtered_read_info[1][1] * 100.0 / stat.n_total
                << '\t' << stat.filtered_read_info[1][1] * 100.0 / stat.n_filtered
                << '\t' << stat.filtered_read_info[1][1] * 100.0 / stat.filtered_read_info[1][4]
                << std::endl;
            filtered_read_info_file << std::fixed << std::setprecision(2) << "high_low_quality_rate"
                << '\t' << stat.filtered_read_info[0][2]
                << '\t' << stat.filtered_read_info[0][2] * 100.0 / stat.n_total
                << '\t' << stat.filtered_read_info[0][2] * 100.0 / stat.n_filtered
                << '\t' << stat.filtered_read_info[0][2] * 100.0 / stat.filtered_read_info[0][4]
                << '\t' << stat.filtered_read_info[1][2] 
                << '\t' << stat.filtered_read_info[1][2] * 100.0 / stat.n_total
                << '\t' << stat.filtered_read_info[1][2] * 100.0 / stat.n_filtered
                << '\t' << stat.filtered_read_info[1][2] * 100.0 / stat.filtered_read_info[1][4]
                << std::endl;
            filtered_read_info_file << std::fixed << std::setprecision(2) << "with_adapter"
                << '\t' << stat.filtered_read_info[0][3]
                << '\t' << stat.filtered_read_info[0][3] * 100.0 / stat.n_total
                << '\t' << stat.filtered_read_info[0][3] * 100.0 / stat.n_filtered
                << '\t' << stat.filtered_read_info[0][3] * 100.0 / stat.filtered_read_info[0][4]
                << '\t' << stat.filtered_read_info[1][3] 
                << '\t' << stat.filtered_read_info[1][3] * 100.0 / stat.n_total
                << '\t' << stat.filtered_read_info[1][3] * 100.0 / stat.n_filtered
                << '\t' << stat.filtered_read_info[1][3] * 100.0 / stat.filtered_read_info[1][4]
                << std::endl;
        }
        for (int i = 0; i < n_end; i++) {
            std::ofstream base_info_file(out_dir.string() + "/Base_distribution_" + std::to_string(i + 1) + ".txt", std::ios_base::out);
            std::ofstream raw_base_quality_info_file(out_dir.string() + "/Base_quality_distribution_raw_" + std::to_string(i + 1) + ".txt", std::ios_base::out);
            std::ofstream clean_base_quality_info_file(out_dir.string() + "/Base_quality_distribution_clean_" + std::to_string(i + 1) + ".txt", std::ios_base::out);

            // base info
            base_info_file << "Len_Pos\tn_raw\t%_raw\tA_raw\t%_raw\tC_raw\t%_raw\tG_raw\t%_raw\tT_raw\t%_raw\tN_raw\t%_raw\tn_raw\t%_raw\tA_clean\t%_clean\tC_clean\t%_clean\tG_clean\t%_clean\tT_clean\t%_clean\tN_clean\t%_clean" << std::endl;
            std::vector<unsigned long> sum_by_base(stat.base_info[i][0].size());
            // unsigned long raw_read_len_sum = std::accumulate(stat.read_len_info[2 * i].begin(), stat.read_len_info[2 * i].end(), 0);
            // unsigned long clean_read_len_sum = std::accumulate(stat.read_len_info[2 * i + 1].begin(), stat.read_len_info[2 * i + 1].end(), 0);
            unsigned long raw_base_sum = 0;
            unsigned long clean_base_sum = 0;
            for (int j = 0; j < stat.base_info[i].size(); j++) {
                base_info_file << std::fixed << std::setprecision(2) << j + 1 << '\t' << stat.read_len_info[2 * i][j] << '\t' << stat.read_len_info[2 * i][j] * 100.0 / stat.n_total;
                for (int k = 0; k < stat.base_info[i][j].size() / 2; k++) {
                    base_info_file << '\t' << std::fixed << std::setprecision(2) << stat.base_info[i][j][k] << '\t' << stat.base_info[i][j][k] * 100.0 / stat.n_total;
                    sum_by_base[k] += stat.base_info[i][j][k];
                    raw_base_sum += stat.base_info[i][j][k];
                }
                base_info_file << std::fixed << std::setprecision(2) << '\t' << stat.read_len_info[2 * i + 1][j] << '\t' << stat.read_len_info[2 * i + 1][j] * 100.0 / stat.n_clean;
                for (int k = stat.base_info[i][j].size() / 2; k < stat.base_info[i][j].size(); k++) {
                    base_info_file << '\t' << std::fixed << std::setprecision(2) << stat.base_info[i][j][k] << '\t' << stat.base_info[i][j][k] * 100.0 / stat.n_clean;
                    sum_by_base[k] += stat.base_info[i][j][k];
                    clean_base_sum += stat.base_info[i][j][k];
                }
                base_info_file << std::endl;
            }
            base_info_file << std::fixed << std::setprecision(2) << "Total\t" << stat.n_total << '\t' << stat.n_total * 100.0 / stat.n_total;
            for (int j = 0; j < sum_by_base.size() / 2; j++) {
                base_info_file << '\t' << std::fixed << std::setprecision(2) << sum_by_base[j] << '\t' << sum_by_base[j] * 100.0 / raw_base_sum;
            }
            base_info_file << std::fixed << std::setprecision(2) << '\t' << stat.n_clean << '\t' << stat.n_clean * 100.0 / stat.n_clean;
            for (int j = sum_by_base.size() / 2; j < sum_by_base.size(); j++) {
                base_info_file << '\t' << std::fixed << std::setprecision(2) << sum_by_base[j] << '\t' << sum_by_base[j] * 100.0 / clean_base_sum;
            }
            base_info_file << std::endl;

            // base quality
            raw_base_quality_info_file << "POS";
            for (int j = 0; j < stat.base_quality_info[2 * i][0].size(); j++) {
                raw_base_quality_info_file << "\tQ" << j;
            }
            raw_base_quality_info_file << "\tMean\t10-quantile\t25-quantile\tMedian\t75-quantile\t90-quantile" << std::endl;
            std::vector<unsigned long> sum_by_quality(stat.base_quality_info[i][0].size());
            for (int j = 0; j < stat.base_quality_info[2 * i].size(); j++) {
                raw_base_quality_info_file << j + 1;
                unsigned long quality_sum = 0;
                unsigned long count = 0;
                int quantile_10 = 0;
                int quantile_25 = 0;
                int median = 0;
                int quantile_75 = 0;
                int quantile_90 = 0;
                for (int k = 0; k < stat.base_quality_info[2 * i][j].size(); k++) {
                    raw_base_quality_info_file << '\t' << stat.base_quality_info[2 * i][j][k];
                    quality_sum += k * stat.base_quality_info[2 * i][j][k];
                    sum_by_quality[k] += stat.base_quality_info[2 * i][j][k];
                    count += stat.base_quality_info[2 * i][j][k];
                    if (quantile_10 == 0 && count > stat.n_total * 0.1) {
                        if (count - stat.base_quality_info[2 * i][j][k] < stat.n_total * 0.1) {
                            quantile_10 = k;
                        }
                        else {
                            quantile_10 = k - 1;
                        }
                    }
                    if (quantile_25 == 0 && count > stat.n_total * 0.25) {
                        if (count - stat.base_quality_info[2 * i][j][k] < stat.n_total * 0.25) {
                            quantile_25 = k;
                        }
                        else {
                            quantile_25 = k - 1;
                        }
                    }
                    if (median == 0 && count > stat.n_total * 0.5) {
                        if (count - stat.base_quality_info[2 * i][j][k] < stat.n_total * 0.5) {
                            median = k;
                        }
                        else {
                            median = k - 1;
                        }
                    }
                    if (quantile_75 == 0 && count > stat.n_total * 0.75) {
                        if (count - stat.base_quality_info[2 * i][j][k] < stat.n_total * 0.75) {
                            quantile_75 = k;
                        }
                        else {
                            quantile_75 = k - 1;
                        }
                    }
                    if (quantile_90 == 0 && count > stat.n_total * 0.9) {
                        if (count - stat.base_quality_info[2 * i][j][k] < stat.n_total * 0.9) {
                            quantile_90 = k;
                        }
                        else {
                            quantile_90 = k - 1;
                        }
                    }
                }
                raw_base_quality_info_file << std::fixed << std::setprecision(2) 
                    << '\t' << quality_sum * 1.0 / stat.n_total 
                    << '\t' << quantile_10
                    << '\t' << quantile_25
                    << '\t' << median
                    << '\t' << quantile_75
                    << '\t' << quantile_90
                    << std::endl;
            }
            raw_base_quality_info_file << "Total";
            for (int j = 0; j < sum_by_quality.size(); j++) {
                raw_base_quality_info_file << '\t' << sum_by_quality[j];
            }
            raw_base_quality_info_file << std::endl;
            raw_base_quality_info_file << "Percent_%";
            unsigned long sum_by_quality_total = std::accumulate(sum_by_quality.begin(), sum_by_quality.end(), 0UL);
            for (int j = 0; j < sum_by_quality.size(); j++) {
                raw_base_quality_info_file << std::fixed << std::setprecision(2) << '\t' << sum_by_quality[j] * 100.0 / sum_by_quality_total;
            }
            raw_base_quality_info_file << std::endl;
            clean_base_quality_info_file << "POS";
            for (int j = 0; j < stat.base_quality_info[2 * i + 1][0].size(); j++) {
                clean_base_quality_info_file << "\tQ" << j;
            }
            clean_base_quality_info_file << "\tMean\t10-quantile\t25-quantile\tMedian\t75-quantile\t90-quantile" << std::endl;
            sum_by_quality = std::vector<unsigned long>(stat.base_quality_info[2 * i + 1][0].size());
            for (int j = 0; j < stat.base_quality_info[2 * i + 1].size(); j++) {
                clean_base_quality_info_file << j + 1;
                unsigned long quality_sum = 0;
                unsigned long count = 0;
                int quantile_10 = 0;
                int quantile_25 = 0;
                int median = 0;
                int quantile_75 = 0;
                int quantile_90 = 0;
                for (int k = 0; k < stat.base_quality_info[2 * i + 1][j].size(); k++) {
                    clean_base_quality_info_file << '\t' << stat.base_quality_info[2 * i + 1][j][k];
                    quality_sum += k * stat.base_quality_info[2 * i + 1][j][k];
                    sum_by_quality[k] += stat.base_quality_info[2 * i + 1][j][k];
                    count += stat.base_quality_info[2 * i + 1][j][k];
                    if (quantile_10 == 0 && count > stat.n_clean * 0.1) {
                        if (count - stat.base_quality_info[2 * i + 1][j][k] < stat.n_clean * 0.1) {
                            quantile_10 = k;
                        }
                        else {
                            quantile_10 = k - 1;
                        }
                    }
                    if (quantile_25 == 0 && count > stat.n_clean * 0.25) {
                        if (count - stat.base_quality_info[2 * i + 1][j][k] < stat.n_clean * 0.25) {
                            quantile_25 = k;
                        }
                        else {
                            quantile_25 = k - 1;
                        }
                    }
                    if (median == 0 && count > stat.n_clean * 0.5) {
                        if (count - stat.base_quality_info[2 * i + 1][j][k] < stat.n_clean * 0.5) {
                            median = k;
                        }
                        else {
                            median = k - 1;
                        }
                    }
                    if (quantile_75 == 0 && count > stat.n_clean * 0.75) {
                        if (count - stat.base_quality_info[2 * i + 1][j][k] < stat.n_clean * 0.75) {
                            quantile_75 = k;
                        }
                        else {
                            quantile_75 = k - 1;
                        }
                    }
                    if (quantile_90 == 0 && count > stat.n_clean * 0.9) {
                        if (count - stat.base_quality_info[2 * i + 1][j][k] < stat.n_clean * 0.9) {
                            quantile_90 = k;
                        }
                        else {
                            quantile_90 = k - 1;
                        }
                    }
                }
                clean_base_quality_info_file << std::fixed << std::setprecision(2) 
                    << '\t' << quality_sum * 1.0 / stat.n_clean 
                    << '\t' << quantile_10
                    << '\t' << quantile_25
                    << '\t' << median
                    << '\t' << quantile_75
                    << '\t' << quantile_90
                    << std::endl;
            }
            clean_base_quality_info_file << "Total";
            for (int j = 0; j < sum_by_quality.size(); j++) {
                clean_base_quality_info_file << '\t' << sum_by_quality[j];
            }
            clean_base_quality_info_file << std::endl;
            clean_base_quality_info_file << "Percent_%";
            sum_by_quality_total = std::accumulate(sum_by_quality.begin(), sum_by_quality.end(), 0UL);
            for (int j = 0; j < sum_by_quality.size(); j++) {
                clean_base_quality_info_file << std::fixed << std::setprecision(2) << '\t' << sum_by_quality[j] * 100.0 / sum_by_quality_total;
            }
            clean_base_quality_info_file << std::endl;

            filtered_read_info_file.close();
            base_info_file.close();
            raw_base_quality_info_file.close();
            clean_base_quality_info_file.close();
        }
    }
}
