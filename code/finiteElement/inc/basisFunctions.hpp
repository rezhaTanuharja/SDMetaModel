/**
 * @file basisFunctions.hpp
 *
 * @brief
 * Contains declarations of one-dimensional basis functions.
 *
 * @author Rezha Adrian Tanuharja <rezhadr@outlook.com>
 * @date 2023-11-02
*/

#ifndef BASIS_FUNCTIONS_DECLARATIONS
#define BASIS_FUNCTIONS_DECLARATIONS

#ifndef ALIAS_DECLARATIONS
#include "alias.hpp"
#endif

#include <algorithm>
#include <array>
#include <exception>


namespace mFEM {


/**
 * @brief
 * A quadratic basis function.
 *
 * @param xi          parametric location in range [-1.0, 1.0]
 * @param funcIndex   index of basis function in (0, 1, or 2)
 * @param diffOrder   differential order
 *
 * @return basis function or its derivative's value at given location
*/
FLOAT quadraticBasis (
  const FLOAT xi, const INT funcIndex, const INT diffOrder
);



} // namespace mFEM


#endif // BASIS_FUNCTIONS_DECLARATIONS

