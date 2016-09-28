#include <unordered_set>

namespace fastq_filter {
    std::string log_title();
    int* check_quality_system(const boost::filesystem::path&);
    float get_base_N_rate(const std::string&);
    float get_average_quality(const std::string&, const int);
    float get_low_quality_rate(const std::string&, const int, const int);
    void load_adapter(boost::filesystem::path&, std::unordered_set<std::string>&);
}
