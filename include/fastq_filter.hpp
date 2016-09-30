#include <unordered_set>

namespace fastq_filter {
    std::string log_title();
    int* check_quality_system(const boost::filesystem::path&);
    float get_base_N_rate(const std::string&);
    float get_average_quality(const std::string&, const int);
    float get_low_quality_rate(const std::string&, const int, const int);
    std::unordered_set<std::string> load_adapter(boost::filesystem::path&);
    void processor(std::vector<boost::filesystem::path>&,
            std::vector<boost::filesystem::path>&,
            std::vector<boost::filesystem::path>&,
            boost::filesystem::path&,
            std::vector< std::unordered_set<std::string> >&,
            int*,
            float*,
            unsigned int*,
            int*);
    void merge(std::vector<boost::filesystem::path>&,
            std::vector<boost::filesystem::path>&,
            boost::filesystem::path&,
            int);
}
