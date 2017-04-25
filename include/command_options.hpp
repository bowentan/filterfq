namespace command_options {
    void print_usage();
    void check_option_dependency(int, const boost::program_options::variables_map&, const char*, ...);
    void check_option_independency(int, const boost::program_options::variables_map&, ...);
    int * parse_trim_param(std::vector<std::string> &, int);
}
