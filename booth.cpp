#include <cmath>
#include <iostream>
#include <vector>

#include "nelder-mead.h"

template <typename T>
T booth(const std::vector<T> &m) {
    T sum = 0;
    auto x = m[0];
    auto y = m[1];
    auto first_term = (x + 2.0 * y - 7.0);
    auto second_term = (2.0 * x + y - 5.0);
    sum += first_term * first_term;
    sum += second_term * second_term;
    return sum;
}

int main() {
    constexpr double tol = 1.0e-6;
    try {
        auto starting_point   = std::vector<double>{0.0, 0.0};
        auto sum = nelder_mead::find_min(booth<double>, starting_point, true, {}, {}, {}, tol, tol);
        std::cout << "Booth solution: ";
        for (auto i : sum)
            std::cout << i << " ";
        std::cout << std::endl;
    } catch (std::exception &e) { std::cout << e.what() << std::endl; }
}
