#include <iostream>
#include <fstream>
#include <unordered_set>
#include <boost/filesystem.hpp>
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
#include <quality_system.hpp>


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

    void load_adapter(boost::filesystem::path& adapter_list, std::unordered_set<std::string>& adapter_read_id) {
        std::ifstream inadapter(adapter_list.c_str(), std::ios_base::in | std::ios_base::binary);
        boost::iostreams::filtering_istream inadapter_decompressor;
        inadapter_decompressor.push(boost::iostreams::gzip_decompressor());
        inadapter_decompressor.push(inadapter);

        std::cout << log_title() << "INFO -- Loading adapter list(s)..." << std::endl;
        std::string line;
        std::vector<std::string> line_splited;
        while (getline(inadapter_decompressor, line)) {
            boost::split(line_splited, line, boost::is_any_of("\t"));
            adapter_read_id.insert(line_splited[0]);
            if (adapter_read_id.size() % 10000 == 0) {
                std::cout << log_title() << "INFO -- "
                    << adapter_read_id.size() << " read ID have been loaded." << std::endl;
            }
        }
        std::cout << log_title() << "INFO -- Totally "
            << adapter_read_id.size() << " have been loaded." << std::endl;
        
        // std::cout << log_title() << "INFO -- Sorting adapter read IDs..." << std::endl;
        // std::sort(adapter_read_id.begin(), adapter_read_id.end());
    }
}
