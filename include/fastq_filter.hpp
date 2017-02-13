#include <boost/filesystem.hpp>
#include <unordered_set>

namespace fastq_filter {
    struct statistic {
        unsigned long n_total;
        unsigned long n_filtered;
        unsigned long n_clean;
        std::vector< std::vector<unsigned long> > read_len_info;                      // 2i: raw, 2i+1: clean, size: 2n_end x max_read_len
        std::vector< std::vector<unsigned long> > filtered_read_info;                 // 0: high N rate, 1: low ave quality, 2: high low-quality rate, 3: adapter, 4: filtered, size: n_end x 5
        std::vector< std::vector< std::vector<unsigned long> > > base_info;           // ACGTN, clean ACGTN, size: n_end x max_read_len x 10
        std::vector< std::vector< std::vector<unsigned long> > > base_quality_info;   // 2i: raw, 2i+1: clean, size: 2n_end x max_read_len x 42

        statistic(int, int);
    };

    std::string log_title();
    int* check_quality_system(const boost::filesystem::path&);
    float get_base_N_rate(const std::string&);
    float get_average_quality(const std::string&, const int);
    float get_low_quality_rate(const std::string&, const int, const int);
    int* get_base_info(const std::string&);
    int* get_base_quality_info(const std::string&, const int);
    std::unordered_set<std::string> load_adapter(boost::filesystem::path&);
    void processor(std::vector<boost::filesystem::path>&,
            std::vector<boost::filesystem::path>&,
            std::vector<boost::filesystem::path>&,
            boost::filesystem::path&,
            std::vector< std::unordered_set<std::string> >&,
            int*,
            float*,
            statistic*,
            int*);
    void merge(std::vector<boost::filesystem::path>&,
            std::vector<boost::filesystem::path>&,
            boost::filesystem::path&,
            int);
    void write_statistic(statistic stat, boost::filesystem::path&);
}
