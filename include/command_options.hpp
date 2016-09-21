namespace command_options {
    void check_option_dependency(int, const boost::program_options::variables_map&, const char*, ...);
    void check_option_independency(int, const boost::program_options::variables_map&, ...);
}
