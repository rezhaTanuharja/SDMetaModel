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


/**
 * @brief
 * Read ASCII file and store in a vector.
 *
 * @param filename  name of file to read
 * @param vector    vector to save data
*/
void loadVector (
  const std::string& filename, std::vector<INT>& vector
);


/**
 * @brief
 * Pad matrix with rows of zeros.
 *
 * @param matrix          matrix to pad
 * @param internalDofIds  IDs of original matrix rows
 * @param boundaryDofIds  IDs of additional zero rows
 *
 * @return the padded matrix
*/
MatrixXF zeroPadRows (
  const MatrixXF& matrix,
  const std::vector<INT>& internalDofIds, 
  const std::vector<INT>& boundaryDofIds
);


/**
 * @brief
 * Combine two vectors using two vectors of IDs.
 *
 * @param Xn              vector 1
 * @param Xb              vector 2
 * @param internalDofIds  IDs of values in Xn
 * @param boundaryDofIds  IDs of values in Xb
 *
 * @return the combined vector
*/
VectorXF combine (
  VectorXF Xn, VectorXF Xb, 
  const std::vector<INT>& internalDofIds, 
  const std::vector<INT>& boundaryDofIds
);


} // namespace mFEM


#endif // UTILITIES_DECLARATIONS

