#include <cmath>
#include <iostream>
#include <vector>

#include "nelder-mead-mine.h"

template <typename T>
T cosine(const std::vector<T> &m) {
    T sum = cos(m[0]);
    return sum;
}

int main() {
    try {
        constexpr double tol = 1.0e-6;
        auto starting_point   = std::vector<double>{0.0};


        auto searcher = nelder_mead_mine::Searcher(
            starting_point, {}, {}, true);
        //searcher.verbose = true;
        searcher.function_tolerance(tol);
        searcher.point_tolerance(tol);
        bool converged{false};
        double sum;
        while(!converged) {
            // request a point
            auto point = searcher.get_next_point();
            // run the function
            sum = cosine(point);
            // report the result
            searcher.report(sum);
            // check for convergence
            converged = searcher.converged();
        }
        std::cout << "Cosine solution: ";
        auto point = searcher.get_res();
        for (auto i : point)
            std::cout << i << " ";
        std::cout << std::endl;
    } catch (std::exception &e) { std::cout << e.what() << std::endl; }
}
