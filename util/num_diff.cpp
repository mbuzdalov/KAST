#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <charconv>

int main(int argc, char *argv[]) {
    if (argc != 3) {
        std::cerr << "Error: expected exactly two arguments, the files to compare" << std::endl;
        return 1;
    }

    std::ifstream l(argv[1]);
    std::ifstream r(argv[2]);

    size_t count_tokens = 0;
    bool any_diff = false;
    std::string ls, rs;
    double lv, rv;
    while (!!(l >> ls) & !!(r >> rs)) {
        ++count_tokens;
        if (ls == rs) {
            continue;
        }
        auto [ptr_l, ec_l] = std::from_chars(ls.data(), ls.data() + ls.size(), lv);
        auto [ptr_r, ec_r] = std::from_chars(rs.data(), rs.data() + rs.size(), rv);
        if (ec_l == std::errc() && ec_r == std::errc()) {
            if (std::abs(lv - rv) > 1e-15) {
                std::cerr << "Token " << count_tokens << ": left = " << lv << ", right = " << rv << ", diff = " << lv - rv << std::endl;
                any_diff = true;
            }
        } else {
            std::cerr << "Token " << count_tokens << ": left = '" << ls << "', right = '" << rs << "'" << std::endl;
            any_diff = true;
        }
    }
    if (l || r) {
        if (l) {
            std::cerr << "Left continues after " << count_tokens << " tokens" << std::endl;
        } else {
            std::cerr << "Right continues after " << count_tokens << " tokens" << std::endl;
        }
        any_diff = true;
    }
    return any_diff ? 1 : 0;
}
