#include <string>

namespace quality_system {
    char zero_quality[5] = {'!', '@', '@', '@', '!'};

    void quality_system_convert(std::string& quality_seq, const int from_sys, const int to_sys) {
        int diff = zero_quality[to_sys] - zero_quality[from_sys];
        for (std::string::iterator c = quality_seq.begin(); c != quality_seq.end(); c++) {
            *c += diff;
            if (*c < zero_quality[to_sys]) {
                *c = zero_quality[to_sys];
            }
            if (to_sys == 3) {
                if (*c < 'B') {
                    *c = 'B';
                }
            }
        }
    }
}
