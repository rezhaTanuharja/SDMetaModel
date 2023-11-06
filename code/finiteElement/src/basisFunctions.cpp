/**
 * @file basisFunctions.cpp
 *
 * @brief
 * Contains implementations of one-dimensional basis functions.
 *
 * @author Rezha Adrian Tanuharja <rezhadr@outlook.com>
 * @date 2023-11-02
*/

#ifndef BASIS_FUNCTIONS_DECLARATIONS
#include "basisFunctions.hpp"
#endif


namespace mFEM {


FLOAT quadraticBasis (
  const FLOAT xi, const INT funcIndex, const INT diffOrder
) {

  #ifdef ADDITIONAL_DEBUG

  if ( xi > 1.0 || xi < -1.0 ) {
    throw std::runtime_error (
      "mFEM::quadraticBasis: parametric location out of range!"
    );
  }

  INT allowableIdx[] = {0, 1, 2};
  INT size = sizeof(allowableIdx) / sizeof(allowableIdx[0]);

  auto it = std::find ( allowableIdx, allowableIdx + size, funcIndex );

  if ( it == allowableIdx + size ) {
    throw std::runtime_error (
      "mFEM::quadraticBasis: invalid function index!"
    );
  }

  if ( diffOrder < 0 ) {
    throw std::runtime_error (
      "mFEM::quadraticBasis: negative differential order!"
    );
  }

  #endif // ADDITIONAL_DEBUG

  if ( diffOrder > 2 ) {
    return 0.0;
  }

  if ( diffOrder > 1 ) {
    FLOAT secondDerivatives[] = {
      1.0, -2.0, 1.0
    };
    return secondDerivatives[funcIndex];
  }

  if ( diffOrder > 0 ) {
    FLOAT firstDerivative[] = {
      xi - 0.5,
      -2.0 * xi,
      xi + 0.5
    };
    return firstDerivative[funcIndex];
  }

  FLOAT basisFunctions[] = {
    0.5 * ( xi * xi - xi ),
    1.0 - ( xi * xi      ),
    0.5 * ( xi * xi + xi )
  };

  return basisFunctions[funcIndex];

} // quadraticBasis


} // namespace mFEM

