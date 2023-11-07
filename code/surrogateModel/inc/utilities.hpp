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

#include <algorithm>
#include <fstream>
#include <random>

#include <Eigen/Dense>


namespace surrogate {


#ifdef HIGH_PRECISION
#define FLOAT long double
#define INT   long int
#else
#define FLOAT double
#define INT   int
#endif // HIGH_PRECISION


using VectorXF = Eigen::Vector<FLOAT,Eigen::Dynamic>;
using MatrixXF = Eigen::Matrix<FLOAT,Eigen::Dynamic,Eigen::Dynamic>;
using DiagMatrixXF = Eigen::DiagonalMatrix<FLOAT,Eigen::Dynamic>;


/**
 * @brief
 * Compute the factorial recursively.
 *
 * @param  n
 * @return n!
*/
INT factorial ( const INT n );


/**
 * @brief
 * Compute the binomial coefficient nCk recursively.
 *
 * @param n the size of set
 * @param k the size of subset
 *
 * @return nCk
*/
INT binomialCoeff ( const INT n, const INT k );


/**
 * @brief
 * Compute the physicist Hermite function recursively.
 *
 * @param index index of Hermite function
 * @param x     argument for the function
 *
 * @return He_{index} ( x )
*/
FLOAT hermiteFunction ( const INT index, const FLOAT x );


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
 * Select a random subset from a given set.
 *
 * @param set   the original set to select from
 * @param size  the size of the random subset
 * @param seed  the seed for the random generator
 *
 * @return a random subset of the original set in a vector
*/
std::vector<INT> randomSubset (
  const std::vector<INT>& set, const INT size, const INT seed
);


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
 * Load vector from a binary file.
 *
 * @param filename  name of file to append to
 * @param vector    vector that stores the read values
*/
void loadVector (
  const std::string& filename, std::vector<FLOAT>& vector
);


/**
 * @brief
 * Load vector from an ASCII file.
 *
 * @param filename  name of file to append to
 * @param vector    vector that stores the read values
*/
void loadVector (
  const std::string& filename, std::vector<INT>& vector
);


} // namespace surrogate


#endif // UTILITIES_DECLARATIONS

