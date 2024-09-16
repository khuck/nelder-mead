#ifndef NELDER_MEAD_CPP_H
#define NELDER_MEAD_CPP_H

#include <cmath>
#include <stdexcept>
#include <vector>
#include <sstream>
#include <cassert>

namespace nelder_mead {

template <typename T> class Vec {
  public:
    Vec() {}
    Vec(unsigned int n) : n(n) {
        val.resize(n, 0);
    }
    Vec(std::initializer_list<T> c) {
        n = c.size();
        val.resize(n);
        std::copy(c.begin().c.end(), val.begin());
    }
    Vec(const Vec &lhs) {
        val = lhs.val;
        n   = lhs.n;
    }
    Vec(const std::vector<T> &lhs) {
        val = lhs;
        n   = lhs.size();
    }
    T operator()(unsigned int idx) const {
        if (idx >= n) {
            throw std::range_error("Element access out of range");
        }
        return val[idx];
    }
    T &operator()(unsigned int idx) {
        if (idx >= n) {
            throw std::range_error("Element access out of range");
        }
        return val[idx];
    }

    Vec operator=(const Vec &rhs) {
        val = rhs.val;
        n   = rhs.n;
        return *this;
    }

    Vec operator=(const std::vector<T> &rhs) {
        val = rhs;
        n   = rhs.size();
        return *this;
    }

    Vec operator+(const Vec &rhs) const {
        Vec lhs(n);
        for (unsigned int i = 0; i < n; i++) {
            lhs.val[i] = val[i] + rhs.val[i];
        }
        return lhs;
    }
    Vec operator-(const Vec &rhs) const {
        Vec lhs(n);
        for (unsigned int i = 0; i < n; i++) {
            lhs.val[i] = val[i] - rhs.val[i];
        }
        return lhs;
    }

    Vec operator/(T rhs) const {
        Vec lhs(n);
        for (unsigned int i = 0; i < n; i++) {
            lhs.val[i] = val[i] / rhs;
        }
        return lhs;
    }
    Vec &operator+=(const Vec &rhs) {
        if (n != rhs.n)
            throw std::invalid_argument(
                "The two vectors must have the same length");
        for (unsigned int i = 0; i < n; i++) {
            val[i] += rhs.val[i];
        }
        return *this;
    }

    Vec &operator-=(const Vec &rhs) {
        if (n != rhs.n)
            throw std::invalid_argument(
                "The two vectors must have the same length");
        for (unsigned int i = 0; i < n; i++) {
            val[i] -= rhs.val[i];
        }
        return *this;
    }

    Vec &operator*=(T rhs) {
        for (unsigned int i = 0; i < n; i++) {
            val[i] *= rhs;
        }
        return *this;
    }

    Vec &operator/=(T rhs) {
        for (unsigned int i = 0; i < n; i++) {
            val[i] /= rhs;
        }
        return *this;
    }
    unsigned int size() const {
        return n;
    }
    unsigned int resize(unsigned int _n) {
        val.resize(_n);
        n = _n;
        return n;
    }
    T length() const {
        T ans = 0;
        for (unsigned int i = 0; i < n; i++) {
            ans += (val[i] * val[i]);
        }
        return std::sqrt(ans);
    }
    std::vector<T> &vec() {
        return val;
    }
    void enforce_min(std::vector<T> limit) {
        // don't exceed the limit!
        for (unsigned int i = 0; i < n; i++) {
            val[i] = std::max(val[i], limit[i]);
        }
    }
    void enforce_max(std::vector<T> limit) {
        // don't exceed the limit!
        for (unsigned int i = 0; i < n; i++) {
            val[i] = std::min(val[i], limit[i]);
        }
    }
    std::string to_string(void) {
        std::stringstream ss;
        ss << "[";
        for (unsigned int i = 0; i < n; i++) {
            ss << val[i];
            ss << (i == n-1 ? "]" : ",");
        }
        return ss.str();
    }
    friend Vec operator*(T a, const Vec &b) {
        Vec c(b.size());
        for (unsigned int i = 0; i < b.size(); i++) {
            c.val[i] = a * b.val[i];
        }
        return c;
    }

  private:
    std::vector<T> val;
    unsigned int   n;
};

template <typename Function, typename T = double>
std::vector<T> find_min(const Function &func,
                        std::vector<T> &initial_point,
                        bool adaptive = false,
                        const std::vector<std::vector<T>> &initial_simplex = {},
                        const std::vector<T> &minimum_limit = {},
                        const std::vector<T> &maximum_limit = {},
                        T tol_fun = 1e-8, T tol_x = 1e-8,
                        unsigned int max_iter      = 1000000,
                        unsigned int max_fun_evals = 100000) {
    unsigned int func_evals_count = 0;
    auto         f                = [&](Vec<T> &p) {
        func_evals_count++;
        return func(p.vec());
    };

    // Getting the dimension of function input
    unsigned int dimension = initial_point.size();
    if (dimension <= 0)
        throw std::invalid_argument(
            "A starting point must have at least one dimension.");

    // Setting parameters
    T alpha, beta, gamma, delta;
    if (adaptive) {
        // Using the results of doi:10.1007/s10589-010-9329-3
        alpha = 1;
        beta  = 1 + 2 / dimension;
        gamma = 0.75 - 1 / (2 * dimension);
        delta = 1 - 1 / dimension;
    } else {
        alpha = 1;
        beta  = 2;
        gamma = 0.5;
        delta = 0.5;
    }
    //std::cout << alpha << " " << beta << " "
        //<< gamma << " " << delta << std::endl;


    // Generate initial simplex
    std::vector<Vec<T>> simplex(dimension + 1);
    if (initial_simplex.empty()) {
        simplex[0] = initial_point;
        for (unsigned int i = 1; i <= dimension; i++) {
            Vec<T> p(initial_point);
            T tau = (p(i - 1) < 1e-6 and p(i - 1) > -1e-6) ? 0.00025 : 0.05;
            p(i - 1) += tau;
            simplex[i] = p;
        }
    } else {
        if (initial_simplex.size() != dimension + 1)
            throw std::invalid_argument(
                "The initial simplex must have dimension + 1 elements");
        for (unsigned int i = 0; i < initial_simplex.size(); i++) {
            simplex[i] = initial_simplex[i];
        }
    }

    std::vector<std::pair<bool, T>> value_cache(dimension + 1);
    for (auto &v : value_cache) {
        v.first = false;
    }
    unsigned int biggest_idx  = 0;
    unsigned int smallest_idx = 0;
    T            biggest_val;
    T            second_biggest_val;
    T            smallest_val;
    auto niter = max_iter;
    while (max_iter--) {
        // Find the points that generate the biggest, second biggest and
        // smallest value
        T val;
        // Does the first point in the simplex have a value? if not, evaluate it
        if (not value_cache[0].first) {
            val                   = f(simplex[0]);
            //std::cout << "first 0 " << simplex[0].to_string() << " = " << val << std::endl;
            value_cache[0].first  = true;
            value_cache[0].second = val;
        } else {
            // if we have a value, get it
            val = value_cache[0].second;
        }
        // the first point in the simplex is our baseline, so set things up accoringly
        biggest_val        = val;
        smallest_val       = val;
        second_biggest_val = val;
        biggest_idx        = 0;
        smallest_idx       = 0;
        // iterate over the other points in the simplex
        for (unsigned int i = 1; i < simplex.size(); i++) {
            T val;
            // does this point in the simplex have a value? if not, evaluate it
            if (not value_cache[i].first) {
                val                   = f(simplex[i]);
                //std::cout << "first " << i << " " << simplex[i].to_string() << " = " << val << std::endl;
                value_cache[i].first  = true;
                value_cache[i].second = val;
            } else {
                // if we have a value, get it.
                val = value_cache[i].second;
            }
            // is this value the biggest?
            if (val > biggest_val) {
                biggest_idx = i;
                biggest_val = val;
            // ...or is it the smallest?
            } else if (val < smallest_val) {
                smallest_idx = i;
                smallest_val = val;
            }
        }

        // at this point, we should know the biggest and smallest value
        // and the simplexes that generated them

        // Calculate the difference of function values and the distance between
        // points in the simplex, so that we can know when to stop the
        // optimization
        T max_val_diff   = 0;
        T max_point_diff = 0;
        for (unsigned int i = 0; i < simplex.size(); i++) {
            T val = value_cache[i].second;
            // find the second biggest value, save it for later reflection
            if (i != biggest_idx and val > second_biggest_val) {
                second_biggest_val = val;
            // is this NOT the smallest value in the current set of simplex points?
            } else if (i != smallest_idx) {
                // how far is this point from the current best, and is it the furthest point from the best?
                if (std::abs(val - smallest_val) > max_val_diff) {
                    max_val_diff = std::abs(val - smallest_val);
                }
                // get the manhattan distance of the vector between this point
                // and the smallest one we have seen so far
                T diff = (simplex[i] - simplex[smallest_idx]).length();
                // how "far" is the point from the smallest point?
                if (diff > max_point_diff) {
                    max_point_diff = diff;
                }
            }
        }
        // have we converged? either by being within tolerance of the best - worst or
        // by being within the tolerance of a point distance?
        if ((max_val_diff <= tol_fun and max_point_diff <= tol_x) or
            (func_evals_count >= max_fun_evals) or (max_iter == 0)) {
            std::vector<T> res = simplex[smallest_idx].vec();
            std::cout << "Converged after " << niter - max_iter << " iterations." << std::endl;
            std::cout << "Total func evaluations: " << func_evals_count << std::endl;
            return res;
        }

        // Calculate the centroid of our current set
        Vec<T> x_bar(dimension);
        for (unsigned int i = 0; i < dimension; i++)
            x_bar(i) = 0;
        for (unsigned int i = 0; i < simplex.size(); i++) {
            if (i != biggest_idx)
                x_bar += simplex[i];
        }
        x_bar /= dimension;

        // Calculate the reflection point
        Vec<T> x_r = x_bar + alpha * (x_bar - simplex[biggest_idx]);
        // enforce limits!
        if (minimum_limit.size() == dimension) x_r.enforce_min(minimum_limit);
        if (minimum_limit.size() == dimension) x_r.enforce_max(maximum_limit);
        // evaluate the reflection point
        T reflection_val = f(x_r);
        //std::cout << "reflection = " << x_r.to_string() << " " << reflection_val << std::endl;
        // is it better than our current best?
        if (reflection_val < smallest_val) {
            // Calculate the Expansion point
            Vec<T> x_e = x_bar + beta * (x_r - x_bar);
            // enforce limits!
            if (minimum_limit.size() == dimension) x_e.enforce_min(minimum_limit);
            if (minimum_limit.size() == dimension) x_e.enforce_max(minimum_limit);
            // evaluate the expansion point
            T expansion_val = f(x_e);
            //std::cout << "expansion = " << x_e.to_string() << " " << expansion_val << std::endl;
            // is the expansion point better than the reflection point?
            if (expansion_val < reflection_val) {
                // replace the worst simplex point with our new expansion point
                simplex[biggest_idx] = x_e;
                value_cache[biggest_idx].second = expansion_val;
            } else {
                // otherwise, replace our worst simplex with our reflection point
                simplex[biggest_idx] = x_r;
                value_cache[biggest_idx].second = reflection_val;
            }
        // is the reflection point worse than our second worst?
        } else if (reflection_val >= second_biggest_val) {
            // Compute the contraction point
            Vec<T> x_c(dimension);
            bool outside = false;
            // is the reflection better than the known worst?
            if (reflection_val < biggest_val) {
                // Outside contraction
                outside = true;
                x_c = x_bar + gamma * (x_r - x_bar);
            // is the reflection worse than the known worst?
            } else if (reflection_val >= biggest_val) {
                // Inside contraction
                x_c = x_bar - gamma * (x_r - x_bar);
            }
            // evaluate the contraction point
            T contraction_val = f(x_c);
            //std::cout << "contraction = " << x_c.to_string() << " " << contraction_val << std::endl;
            // is the contraction better than the reflection or the known worst?
            if ((outside and contraction_val <= reflection_val) or
                (not outside and contraction_val <= biggest_val)) {
                // replace the known worst with the contraction
                simplex[biggest_idx] = x_c;
                value_cache[biggest_idx].second = contraction_val;
            } else {
                // Shrinking
                //std::cout << "shrinking" << std::endl;
                // we take the whole simplex, and move every point towards the current best candidate
                for (unsigned int i = 0; i < dimension; i++) {
                    if (i != smallest_idx) {
                        simplex[i] =
                            simplex[smallest_idx] +
                            delta * (simplex[i] - simplex[smallest_idx]);
                        value_cache[i].first = false;
                    }
                }
            }
        } else {
            // Reflection is good enough
            // replace the worst point with the reflection point.
            simplex[biggest_idx] = x_r;
            value_cache[biggest_idx].second = reflection_val;
        }
    }
    // if we hit max iterations, just return the best we found so far
    std::vector<T> res = simplex[smallest_idx].vec();
    return res;
}

}; // namespace nelder_mead
#endif
