#include <cmath>
#include <iostream>
#include <vector>

#include "nelder-mead.h"

template <typename T>
T squared(const std::vector<T> &m) {
    T sum = 0;
    for (unsigned int i = 0; i < m.size(); i++) {
        sum += m[i] * m[i];
    }
    return sum;
}

int main() {
    try {
        double tol = 1.0e-6;
        auto starting_point   = std::vector<double>{1.0};
        auto minimum_limit   = std::vector<double>{-2.0};
        auto maximum_limit   = std::vector<double>{2.0};
        auto searcher = nelder_mead::Searcher<double>(
            starting_point, minimum_limit, maximum_limit);
        searcher.function_tolerance(tol);
        searcher.point_tolerance(tol);
        double sum;
        // check for convergence
        while(!searcher.converged()) {
            // request a point
            auto point = searcher.get_next_point();
            // run the function
            sum = squared(point);
            // report the result
            searcher.report(sum);
        }
        std::cout << "Squared solution: ";
        auto point = searcher.get_res();
        for (auto i : point)
            std::cout << i << " ";
        std::cout << std::endl;
    } catch (std::exception &e) { std::cout << e.what() << std::endl; }
}
