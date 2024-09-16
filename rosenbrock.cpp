#include <cmath>
#include <iostream>
#include <vector>

#include "nelder-mead.h"

// The rosenbrock function
// https://en.wikipedia.org/wiki/Rosenbrock_function
template <typename T>
T rosen(const std::vector<T> &m) {
    T res = 0;
    for (unsigned int i = 0; i < m.size() - 1; i++) {
        res += 100 * (m[i + 1] - m[i] * m[i]) * (m[i + 1] - m[i] * m[i]) +
               (1 - m[i]) * (1 - m[i]);
    }
    return res;
}

int main() {
    try {

        auto starting_point   = std::vector<double>{1, 2, 3, 4, 5};
        auto searcher = nelder_mead_mine::Searcher<double>(
            starting_point);
        double sum;
        // check for convergence
        while(!searcher.converged()) {
            // request a point
            auto point = searcher.get_next_point();
            // run the function
            sum = rosen(point);
            // report the result
            searcher.report(sum);
        }
        std::cout << "Rosen solution: ";
        auto point = searcher.get_res();
        for (auto i : point)
            std::cout << i << " ";
        std::cout << std::endl;
    } catch (std::exception &e) { std::cout << e.what() << std::endl; }
}
