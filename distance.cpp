#include <cmath>
#include <iostream>
#include <vector>

#include "nelder-mead.h"

template <typename T>
T distance(const std::vector<T> &m) {
    T sum = 0;
    for (unsigned int i = 0; i < m.size(); i++) {
        auto tmp = std::abs(m[i]) - i;
        sum += std::abs(tmp * tmp);
    }
    return sum;
}


int main() {
    try {
        auto starting_point   = std::vector<double>{4,4,4,4,4,4,4};
        auto minimum_limit   = std::vector<double>{0,0,0,0,0,0,0};
        auto maximum_limit   = std::vector<double>{7,7,7,7,7,7,7};
        std::array<double,2> tol{0.1,0.01};
        auto searcher = nelder_mead::Searcher<double>(
            starting_point, minimum_limit, maximum_limit, false);
            //starting_point);
        searcher.function_tolerance(tol[0]);
        searcher.point_tolerance(tol[1]);
        bool converged{false};
        while(!converged) {
            // request a point
            auto point = searcher.get_next_point();
            // run the function
            auto sum = distance(point);
            // report the result
            searcher.report(sum);
            // check for convergence
            converged = searcher.converged();
        }
        std::cout << "Distance solution: ";
        auto point = searcher.get_res();
        for (auto i : point)
            std::cout << std::round(i) << " ";
        std::cout << std::endl;
    } catch (std::exception &e) { std::cout << e.what() << std::endl; }
}
