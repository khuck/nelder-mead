#include <cmath>
#include <iostream>
#include <vector>

#include "nelder-mead.h"

template <typename T>
T beale(const std::vector<T> &m) {
    T sum = 0;
    auto x = m[0];
    auto y = m[1];
    auto first_term = (1.5 - x + x * y);
    auto second_term = (2.25 - x + x * y * y);
    auto third_term = (2.625 - x + x * y * y * y);
    sum += first_term * first_term;
    sum += second_term * second_term;
    sum += third_term * third_term;
    return sum;
}

int main() {
    constexpr double tol = 1.0e-6;
    try {
        auto starting_point   = std::vector<double>{0.0, 0.0};
        auto sum = nelder_mead::find_min(beale<double>, starting_point, true, {}, {}, {}, tol, tol);
        std::cout << "Beale solution: ";
        for (auto i : sum)
            std::cout << i << " ";
        std::cout << std::endl;
    } catch (std::exception &e) { std::cout << e.what() << std::endl; }
}
