#include <boost/program_options.hpp>

namespace command_options {
    void inter_dependency(const boost::program_options::variables_map& vm, const char* for_what, const char* required_option) {
        if (!vm.count("help")) {
            if (vm.count(for_what) && !vm[for_what].defaulted())
                if (!vm.count(required_option) && !vm[required_option].defaulted()) {
                    throw std::logic_error(std::string("Option '") + for_what + "' requires option '" + required_option + "'.");
                }
        }
    }

    void inter_independency(const boost::program_options::variables_map& vm, const char* opt1, const char* opt2) {
        if (!vm.count("help")) {
            if (!vm[opt1].defaulted() && !vm[opt2].defaulted()) {
                if (vm.count(opt1) && vm.count(opt2)) {
                    throw std::logic_error(std::string("Conflicting options '") + opt1 + "' and '" + opt2 + "'.");
                }

                if (!vm.count(opt1) && !vm.count(opt2)) {
                    throw std::logic_error(std::string("One of options '") + opt1 + "' and '" + opt2 + "' is requred but missing.");
                }
            }
        }
    }
}
