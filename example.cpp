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

template <typename T>
T squared(const std::vector<T> &m) {
    T sum = 0;
    for (unsigned int i = 0; i < m.size(); i++) {
        sum += m[i] * m[i];
    }
    //std::cout << m[0]  << " = " << sum << std::endl;
    return sum;
}

template <typename T>
T cosine(const std::vector<T> &m) {
    T sum = cos(m[0]);
    //std::cout << m[0]  << " = " << sum << std::endl;
    return sum;
}

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
    //std::cout << m[0]  << " = " << sum << std::endl;
    return sum;
}

template <typename T>
T booth(const std::vector<T> &m) {
    T sum = 0;
    auto x = m[0];
    auto y = m[1];
    auto first_term = (x + 2.0 * y - 7.0);
    auto second_term = (2.0 * x + y - 5.0);
    sum += first_term * first_term;
    sum += second_term * second_term;
    //std::cout << "[" << m[0]  << "," << m[1] << "] = " << sum << std::endl;
    return sum;
}

template <typename T>
T distance(const std::vector<T> &m) {
    T sum = 0;
    for (unsigned int i = 0; i < m.size(); i++) {
        auto tmp = std::abs(m[i]) - i;
        sum += std::abs(tmp * tmp);
    }
    //std::cout << m[0]  << " = " << sum << std::endl;
    return sum;
}


int main() {
    try {
#if 1
        {
        auto starting_point   = std::vector<double>{1, 2, 3, 4, 5};
        auto starting_simplex = std::vector<std::vector<double>>{
            {2, 3, 4, 5, 1}, {3, 4, 5, 1, 2}, {4, 5, 1, 2, 3},
            {5, 1, 2, 3, 4}, {1, 2, 3, 4, 5}, {0, 0, 0, 0, 0}};
        auto res = nelder_mead::find_min(rosen<double>, starting_point, true,
                                         starting_simplex);
        std::cout << "Rosen solution: ";
        for (auto i : res)
            std::cout << i << " ";
        std::cout << std::endl;
        }

        {
        auto starting_point   = std::vector<double>{7,6,5,4,3,2,1};
        std::array<double,2> tol{0.1,0.01};
        auto sum = nelder_mead::find_min(distance<double>, starting_point, true,
            {}, {}, {},
            tol[0],
            tol[1]);
        std::cout << "Distance solution: ";
        for (auto i : sum)
            std::cout << int(std::round(i)) << " ";
        std::cout << std::endl;
        }

        double tol = 1.0e-6;

        {
        auto starting_point   = std::vector<double>{1.0};
        auto sum = nelder_mead::find_min(squared<double>, starting_point, true, {}, {}, {}, tol, tol);
        std::cout << "Squared solution: ";
        for (auto i : sum)
            std::cout << i << " ";
        std::cout << std::endl;
        }
#endif
        {
        double tol = 1.0e-6;
        auto starting_point   = std::vector<double>{0.0};
        auto sum = nelder_mead::find_min(cosine<double>, starting_point, true, {}, {}, {}, tol, tol);
        std::cout << "Cosine solution: ";
        for (auto i : sum)
            std::cout << i << " ";
        std::cout << std::endl;
        }
#if 1
        {
        auto starting_point   = std::vector<double>{0.0, 0.0};
        auto sum = nelder_mead::find_min(beale<double>, starting_point, true, {}, {}, {}, tol, tol);
        std::cout << "Beale solution: ";
        for (auto i : sum)
            std::cout << i << " ";
        std::cout << std::endl;
        }

        {
        auto starting_point   = std::vector<double>{0.0, 0.0};
        auto sum = nelder_mead::find_min(booth<double>, starting_point, true, {}, {}, {}, tol, tol);
        std::cout << "Booth solution: ";
        for (auto i : sum)
            std::cout << i << " ";
        std::cout << std::endl;
        }
#endif

    } catch (std::exception &e) { std::cout << e.what() << std::endl; }
}
