/**
 * @file utilities.hpp
 *
 * @brief
 * Contains declarations of various useful functions.
 *
 * @author Rezha Adrian Tanuharja <rezhadr@outlook.com>
 * @date 2023-11-02
*/

#ifndef UTILITIES_DECLARATIONS
#define UTILITIES_DECLARATIONS

#ifndef ALIAS_DECLARATIONS
#include "alias.hpp"
#endif

#include <array>
#include <cmath>
#include <exception>
#include <fstream>
#include <random>
#include <string>
#include <vector>


namespace mFEM {


/**
 * @brief
 * Gives Gauss integration points and weights.
 *
 * @param n number of integrations points
 *
 * @return vector of integration points and weights
*/
IntegrationPoints gaussQuadrature ( const INT n );


/**
 * @brief
 * Sample a random variable from a normal distribution.
 *
 * @param mean    mean of the underlying normal distribution
 * @param stdDev  std deviation of the underlying normal distribution
 * @param seed    seed for the random generator
 *
 * @return a realization of random variable
*/
FLOAT normalVariable ( const FLOAT mean, const FLOAT stdDev, const INT seed );


/**
 * @brief
 * Append vector to a binary file.
 *
 * @param filename  name of file to append to
 * @param vector    vector of FLOAT to append
*/
void saveVector (
  const std::string& filename, const std::vector<FLOAT>& vector
);


/**
 * @brief
 * Append vector to an ASCII file.
 *
 * @param filename  name of file to append to
 * @param vector    vector of INT to append
*/
void saveVector (
  const std::string& filename, const std::vector<INT>& vector
);


} // namespace mFEM


#endif // UTILITIES_DECLARATIONS

