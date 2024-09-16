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
    try {
        constexpr double tol = 1.0e-6;
        auto starting_point   = std::vector<double>{0.0, 0.0};
        auto searcher = nelder_mead::Searcher<double>(
            starting_point);
        searcher.function_tolerance(tol);
        searcher.point_tolerance(tol);
        double sum;
        // check for convergence
        while(!searcher.converged()) {
            // request a point
            auto point = searcher.get_next_point();
            // run the function
            sum = booth(point);
            // report the result
            searcher.report(sum);
        }
        std::cout << "Booth solution: ";
        auto point = searcher.get_res();
        for (auto i : point)
            std::cout << i << " ";
        std::cout << std::endl;
    } catch (std::exception &e) { std::cout << e.what() << std::endl; }
}
