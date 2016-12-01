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

    std::string log_title() {return "[filterfq | " + to_simple_string(boost::posix_time::second_clock::local_time()) + "] ";}

    int* check_quality_system(const boost::filesystem::path& filepath) {
        int* results = new int[3];
        std::ifstream infile (filepath.string(), std::ios_base::in | std::ios_base::binary);
        boost::iostreams::filtering_istream decompressor;
        decompressor.push(boost::iostreams::gzip_decompressor());
        decompressor.push(infile);
        std::string line;

        unsigned char min = '~';
        unsigned char max = '!';
        int n_read = 1;
        std::pair<std::string::iterator, std::string::iterator> extremum;
        boost::random::mt19937 gen;
        boost::random::uniform_int_distribution<> dist(1, 10);
        while (n_read <= 100000 && getline(decompressor, line)) {
            for (int i = 0; i < 3; i++)
                getline(decompressor, line);
            if (dist(gen) >= 5) 
                continue;

            extremum = boost::minmax_element(line.begin(), line.end());
            // std::cout << line << std::endl;
            // std::cout << *extremum.first << std::endl;
            // std::sort(line.begin(), line.end());
            // std::cout << line << std::endl;
            // std::cout << *extremum.first << std::endl;
            // std::cout << *extremum.first << " " << *extremum.second << std::endl;
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
        // std::cout << "min quality is \"" << min << "\" and max quality is \"" << max << "\"" << std::endl;
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
        return results;
    }
    
    float get_base_N_rate(const std::string& read_seq) {
        return std::count(read_seq.begin(), read_seq.end(), 'N') * 1.0 / read_seq.length();
    }

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
        
        // std::cout << log_title() << "INFO -- Sorting adapter read IDs..." << std::endl;
        // std::sort(adapter_read_id_list.begin(), adapter_read_id_list.end());
    }

    void processor(std::vector<boost::filesystem::path>& infiles,
            std::vector<boost::filesystem::path>& clean_outfiles,
            std::vector<boost::filesystem::path>& dropped_outfiles,
            boost::filesystem::path& tmp_dir,
            std::vector< std::unordered_set<std::string> >& adapter_read_id_lists,
            int* param_int,
            float* param_float,
            unsigned int* counter,
            int* thread_info) {
        int stepsize = thread_info[0];
        int n_thread = thread_info[1];
        int thread = thread_info[2];
        
        int min_base_quality = param_int[0];
        int raw_quality_sys = param_int[1];
        int clean_quality_sys = param_int[2];
        float max_base_N_rate = param_float[0];
        float min_ave_quality = param_float[1];
        float max_low_quality_rate = param_float[2];

        if (infiles.size() == 1) {
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
            
            unsigned int local_counter = 0;
            std::string read_id_line, read_line, plus_line, quality_line;
            for (int i = 0; i < stepsize * thread; i++) {
                for (int j = 0; j < 4; j++) {
                    getline(infq_decompressor, read_id_line);
                }
            }

            bool is_filtered;
            while (getline(infq_decompressor, read_id_line)) {
                getline(infq_decompressor, read_line);
                getline(infq_decompressor, plus_line);
                getline(infq_decompressor, quality_line);

                is_filtered = false;
                // cout << read_id_line << endl
                //     << read_line << endl 
                //     << plus_line << endl
                //     << quality_line << endl;
                // cout << get_base_N_rate(read_line) << " "
                //     << get_average_quality(quality_line, raw_quality_sys) << " "
                //     << get_low_quality_rate(quality_line, raw_quality_sys, min_base_quality) << endl;
                mutex.lock();
                if (adapter_read_id_lists.size() != 0) {
                    std::string s = read_id_line.substr(1, read_id_line.size() - 1);
                    if (adapter_read_id_lists[0].count(s) > 0) {
                        counter[5]++;
                        if (!is_filtered) {
                            counter[1]++;
                            is_filtered = true;
                        }
                    }
                }
                if (get_base_N_rate(read_line) > max_base_N_rate) {
                    counter[2]++;
                    if (!is_filtered) {
                        counter[1]++;
                        is_filtered = true;
                    }
                }
                if (get_average_quality(quality_line, raw_quality_sys) < min_ave_quality) {
                    counter[3]++;
                    if (!is_filtered) {
                        counter[1]++;
                        is_filtered = true;
                    }
                }
                if (get_low_quality_rate(quality_line, raw_quality_sys, min_base_quality) > max_low_quality_rate) {
                    counter[4]++;
                    if (!is_filtered) {
                        counter[1]++;
                        is_filtered = true;
                    }
                }

                counter[0]++;
                if (counter[0] % 50000 == 0) {
                    std::cout << log_title() << "INFO " 
                        << std::fixed
                        << std::setw(12) << std::setprecision(6) << counter[0] << " | "
                        << std::setw(12) << std::setprecision(6) << counter[1] << " | "
                        << std::setw(12) << std::setprecision(6) << counter[1] * 100.0 / counter[0] << " | "
                        << std::setw(12) << std::setprecision(6) << counter[2] << " | "
                        << std::setw(12) << std::setprecision(6) << counter[3] << " | "
                        << std::setw(12) << std::setprecision(6) << counter[4] << " | "
                        << std::setw(12) << std::setprecision(6) << counter[5] << " | "
                        << std::endl;
                }
                mutex.unlock();

                if (!is_filtered) {
                    if (raw_quality_sys != clean_quality_sys) {
                        quality_system::quality_system_convert(quality_line, raw_quality_sys, clean_quality_sys);
                    }
                    clean_outfq_compressor << read_id_line << std::endl
                        << read_line << std::endl
                        << plus_line << std::endl
                        << quality_line << std::endl;
                }
                else {
                    dropped_outfq_compressor << read_id_line << std::endl
                        << read_line << std::endl
                        << plus_line << std::endl
                        << quality_line << std::endl;
                }

                local_counter++;
                if (local_counter % stepsize == 0) {
                    for (int i = 0; i < stepsize * (n_thread - 1); i++) {
                        for (int j = 0; j < 4; j++) {
                            getline(infq_decompressor, read_id_line);
                        }
                    }
                }
            }
            close(infq_decompressor, std::ios_base::in);
            close(clean_outfq_compressor, std::ios_base::out);
            close(dropped_outfq_compressor, std::ios_base::out);
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
            
            unsigned int local_counter = 0;
            std::string read_id_line1, read_id_line2,
                read_line1, read_line2, 
                plus_line1, plus_line2,
                quality_line1, quality_line2;
            for (int i = 0; i < stepsize * thread; i++) {
                for (int j = 0; j < 4; j++) {
                    getline(infq1_decompressor, read_id_line1);
                    getline(infq2_decompressor, read_id_line2);
                }
            }

            bool is_filtered;
            while (getline(infq1_decompressor, read_id_line1) && getline(infq2_decompressor, read_id_line2)) {
                getline(infq1_decompressor, read_line1);
                getline(infq2_decompressor, read_line2);
                getline(infq1_decompressor, plus_line1);
                getline(infq2_decompressor, plus_line2);
                getline(infq1_decompressor, quality_line1);
                getline(infq2_decompressor, quality_line2);

                is_filtered = false;
                // cout << read_id_line << endl
                //     << read_line << endl 
                //     << plus_line << endl
                //     << quality_line << endl;
                // cout << get_base_N_rate(read_line) << " "
                //     << get_average_quality(quality_line, raw_quality_sys) << " "
                //     << get_low_quality_rate(quality_line, raw_quality_sys, min_base_quality) << endl;
                mutex.lock();
                if (adapter_read_id_lists.size() != 0) {
                    std::string s1 = read_id_line1.substr(1, read_id_line1.size() - 1);
                    std::string s2 = read_id_line2.substr(1, read_id_line2.size() - 1);
                    if (adapter_read_id_lists[0].count(s1) > 0 || adapter_read_id_lists[1].count(s2) > 0) {
                        counter[5]++;
                        if (!is_filtered) {
                            counter[1]++;
                            is_filtered = true;
                        }
                    }
                }
                if (get_base_N_rate(read_line1) > max_base_N_rate || get_base_N_rate(read_line2) > max_base_N_rate) {
                    counter[2]++;
                    if (!is_filtered) {
                        counter[1]++;
                        is_filtered = true;
                    }
                }
                if (get_average_quality(quality_line1, raw_quality_sys) < min_ave_quality || get_average_quality(quality_line2, raw_quality_sys) < min_ave_quality) {
                    counter[3]++;
                    if (!is_filtered) {
                        counter[1]++;
                        is_filtered = true;
                    }
                }
                if (get_low_quality_rate(quality_line1, raw_quality_sys, min_base_quality) > max_low_quality_rate || get_low_quality_rate(quality_line2, raw_quality_sys, min_base_quality) > max_low_quality_rate) {
                    counter[4]++;
                    if (!is_filtered) {
                        counter[1]++;
                        is_filtered = true;
                    }
                }

                counter[0]++;
                if (counter[0] % 50000 == 0) {
                    std::cout << log_title() << "INFO " 
                        << std::fixed
                        << std::setw(12) << std::setprecision(6) << counter[0] << " | "
                        << std::setw(12) << std::setprecision(6) << counter[1] << " | "
                        << std::setw(12) << std::setprecision(6) << counter[1] * 100.0 / counter[0] << " | "
                        << std::setw(12) << std::setprecision(6) << counter[2] << " | "
                        << std::setw(12) << std::setprecision(6) << counter[3] << " | "
                        << std::setw(12) << std::setprecision(6) << counter[4] << " | "
                        << std::setw(12) << std::setprecision(6) << counter[5] << " | "
                        << std::endl;
                }
                mutex.unlock();

                if (!is_filtered) {
                    if (raw_quality_sys != clean_quality_sys) {
                        quality_system::quality_system_convert(quality_line1, raw_quality_sys, clean_quality_sys);
                        quality_system::quality_system_convert(quality_line2, raw_quality_sys, clean_quality_sys);
                    }
                    clean_outfq1_compressor << read_id_line1 << std::endl
                        << read_line1 << std::endl
                        << plus_line1 << std::endl
                        << quality_line1 << std::endl;
                    clean_outfq2_compressor << read_id_line2 << std::endl
                        << read_line2 << std::endl
                        << plus_line2 << std::endl
                        << quality_line2 << std::endl;
                }
                else {
                    dropped_outfq1_compressor << read_id_line1 << std::endl
                        << read_line1 << std::endl
                        << plus_line1 << std::endl
                        << quality_line1 << std::endl;
                    dropped_outfq2_compressor << read_id_line2 << std::endl
                        << read_line2 << std::endl
                        << plus_line2 << std::endl
                        << quality_line2 << std::endl;
                }

                local_counter++;
                if (local_counter % stepsize == 0) {
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
        }
    }

    void merge_single(boost::filesystem::path& outfile, boost::filesystem::path& tmp_dir, int n_thread) {
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

    // void merge(std::vector<boost::filesystem::path>& clean_outfiles,
    //         std::vector<boost::filesystem::path>& dropped_outfiles,
    //         boost::filesystem::path& tmp_dir,
    //         int n_thread) {
    //     std::ofstream clean_out[clean_outfiles.size()];
    //     std::ofstream dropped_out[dropped_outfiles.size()];
    //     for (int i = 0; i < clean_outfiles.size(); i++) {
    //         clean_out[i].open(clean_outfiles[i].string(), std::ios_base::out | std::ios_base::binary | std::ios_base::app);
    //         dropped_out[i].open(dropped_outfiles[i].string(), std::ios_base::out | std::ios_base::binary | std::ios_base::app);
    //     }
    //     
    //     for (int i = 0; i < clean_outfiles.size(); i++) {
    //         for (int j = 0; j < n_thread; j++) {
    //             std::string clean_filename = (tmp_dir / clean_outfiles[i].filename()).string() + "." + std::to_string(j) + ".tmp";
    //             std::string dropped_filename = (tmp_dir / dropped_outfiles[i].filename()).string() + "." + std::to_string(j) + ".tmp";
    //             std::ifstream clean_file(clean_filename, std::ios_base::in | std::ios_base::binary);
    //             std::ifstream dropped_file(dropped_filename, std::ios_base::in | std::ios_base::binary);

    //             clean_out[i] << clean_file.rdbuf();
    //             dropped_out[i] << dropped_file.rdbuf();

    //             clean_file.close();
    //             dropped_file.close();

    //             boost::filesystem::remove(clean_filename);
    //             boost::filesystem::remove(dropped_filename);
    //         }
    //     }
    //     
    //     for (int i = 0; i < clean_outfiles.size(); i++) {
    //         clean_out[i].close();
    //         dropped_out[i].close();
    //     }
    // }
}
