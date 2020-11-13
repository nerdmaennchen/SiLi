#pragma once

/*! Overview
 * \page
 *
 * SiLi is a great autostorage only library (no heap allocation).
 *
 * It uses **c++20** and its purpose is to be used on embedded system
 *
 * The main idea is to use fixed size matrices and fixed size
 * views onto these matrices.
 *
 * \creategroup 1 Pages
 * \creategroup 2 Concepts
 * \creategroup 3 Classes
 * \creategroup 4 Free Matrix Function
 * \creategroup 4 Matrix Operations
 * \creategroup 5 Free Vector Function
 */


/*! Examples
 * \page
 * Example of Matrix Multiplication:
 * \code
 *   #include <SiLi/SiLi.h>
 *   #include <SiLi/ostream.h>
 *
 *   int main() {
 *   auto a = SiLi::Matrix{{{3, 4, 7},
 *                          {5, 6, 8}}};
 *   auto b = SiLi::Matrix{{{1, 2},
 *                          {4, 5},
 *                          {3, 6}};
 *   auto c = a * b; // c is of type SiLi::Matrix<2, 2, int>
 *   std::cout << c << "\n"; // prints {{40, 68},
 *                           //         {53, 88}}
 *   }
 * \endcode
 */


#include "Matrix.h"
#include "View.h"
#include "operations.h"
#include "Iterator.h"
