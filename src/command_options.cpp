#include <boost/program_options.hpp>
#include <cstdarg>

namespace command_options {
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
}
