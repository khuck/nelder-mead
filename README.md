# Modified Nelder Mead in C++

This implementation of [Nelder Mead
method](https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method) was inspired
by the implementaiton at https://github.com/YibaiMeng/nelder-mead, however has
been modified for use as a two-step search strategy. The modification allows us
to request a new point in the simplex, use those settings elsewhere in a
program, and then report the resulting output. This removes the main `while`
loop in the previous implementation, and replaces it with a simple state
machine that allows us to iterate through the stages of SIMPLEX, REFLECTION,
EXPANSION, and CONTRACTION as necessary.

I am grateful to the article ["Breaking down the Nelder Mead
algorithm"](https://brandewinder.com/2022/03/31/breaking-down-Nelder-Mead/) for
a great explanation of the algorithm with helpful diagrams.

The original is a simple C++ implementation of Nelder Mead, a numerical method
to find the minimum or maximum of an objective function in a multidimensional
space. That implementation uses the Nelder Mead method in [Implementing the
Nelder-Mead simplex algorithm with adaptive
parameters](https://link.springer.com/article/10.1007/s10589-010-9329-3) by Gao
and Han, which makes a modification to improve convergence in higher
dimensions. The original implementation has been retained in nelder-mead-old.h
as a reference.

This is a header only library. To use this library, simply include the
`nelder_mead.h` file in your program. The library defines one public class,
`Searcher`, with four main methods: a constructor, `get_new_point`, `report`,
and `converged`.

```cpp
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
```
The output is
```
1 3 0.75 0
Converged after 59 iterations.
Total func evaluations: 118
Cosine solution: 3.14159
```
See [`example.cpp`](example.cpp) for examples using the old implementation, and the following test programs to see the new implementation:
* [`beale.cpp`](beale.cpp)
* [`booth.cpp`](booth.cpp)
* [`cosine.cpp`](cosine.cpp)
* [`distance.cpp`](distance.cpp)
* [`rosenbrock.cpp`](rosenbrock.cpp)
* [`squared.cpp`](squared.cpp)

See the paper for a detailed description of the algorithm, including the termination criteria.

This project is released into the Public Domain, as was the original implementation.
