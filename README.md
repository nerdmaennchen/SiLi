# SiLi
A simple and basic linear algebra library for c++ using little code.

Full Reference [here](https://sgssgene.github.io/SiLi.html)

# Features
* Compile time matrices
* no heap allocations
* exchangeable datatype
* Matrix operations:
  * Matrix operations: multiplication, addition, subtraction, negation, assignment
  * Element wise operations: multiplication, assignment
  * views on matrices
  * determinant
  * inverse()
  * norm()
  * transpose (as a view)
  * diagonal access (as a view)
  * iteration over elements, rows or columns possible
  * join_rows()/join_cols()
  * abs()
  * sum()
  * cross() (for 3x1 Matrices)


# Example
```
#include <SiLi/SiLi.h>
#include <SiLi/ostream.h>

int main() {
    // m1 is of type Matrix<3, 2, double>
    auto m1 = SiLi::Matrix{{{1., 0.},
                            {0., 1.},
                            {3., 4.}}};

    // m2 is of type Matrix<2, 3, double>
    auto m2 = SiLi::Matrix{{{1., 0., 3.},
                            {0., 1., 4.}}};

    // m3 is of type Matrix<3, 3, double>
    auto m3 = m1 * m2;

    // m4 is of type Matrix<2, 2, double>
    auto m4 = m2 * m1;

    // print matrix
    std::cout << "m1:\n" << m1 << "\n\n";
    std::cout << "m2:\n" << m2 << "\n\n";
    std::cout << "m3:\n" << m3 << "\n\n";
    std::cout << "m4:\n" << m4 << "\n\n";

    // print matrix element by element
    for (int row = 0; row < rows(m1); ++row) {
        for (int col = 0; col < cols(m1); ++col) {
            std::cout << m1(row, col) << "\n";
        }
    }
}
```

