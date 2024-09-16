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
            sum = beale(point);
            // report the result
            searcher.report(sum);
        }
        std::cout << "Beale solution: ";
        auto point = searcher.get_res();
        for (auto i : point)
            std::cout << i << " ";
        std::cout << std::endl;
    } catch (std::exception &e) { std::cout << e.what() << std::endl; }
}
