namespace fastq_filter {
    int* check_quality_system(const boost::filesystem::path&);
    float get_base_N_rate(const std::string&);
    float get_average_quality(const std::string&, const int);
    float get_low_quality_rate(const std::string&, const int, const int);
}
