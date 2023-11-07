/**
 * @file basisFunctions.hpp
 *
 * @brief
 * Contains declarations of basis functions for PCE and RPCE.
 *
 * @author Rezha Adrian Tanuharja <rezhadr@outlook.com>
 * @date 2023-11-02
*/

#ifndef BASIS_FUNCTIONS_DECLARATIONS
#define BASIS_FUNCTIONS_DECLARATIONS

#ifndef UTILITIES_DECLARATIONS
#include "utilities.hpp"
#endif


namespace surrogate {


/**
 * @brief
 * Generate sets of indices of functions.
 * Sets are sorted from lowest sum of indices.
 *
 * @param dimension number of indices in each set
 * @param maxSum    the largest sum of indices in a set
 *
 * @return vector of indices
*/
std::vector<INT> generateIndices (
  const INT dimension, const INT maxSum
);


/**
 * @brief
 * Generate basis functions from products of probabilist Hermite functions.
 *
 * @param indices   indices of Hermite functions in sequence
 * @param args      arguments for the Hermite functions
 * @param dimension number of Hermite functions to multiply for a basis
 *
 * @return vector of Hermite basis functions
*/
std::vector<FLOAT> hermiteBasis (
  const std::vector<INT>& indices, 
  const std::vector<FLOAT>& args, 
  const INT dimension
);


} // namespace surrogate


#endif // BASIS_FUNCTIONS_DECLARATIONS

