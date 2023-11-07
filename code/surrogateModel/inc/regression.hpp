/**
 * @file regression.hpp
 *
 * @brief
 * Contains declarations of various functions for regression.
 *
 * @author Rezha Adrian Tanuharja <rezhadr@outlook.com>
 * @date 2023-11-02
*/

#ifndef REGRESSION_DECLARATIONS
#define REGRESSION_DECLARATIONS

#ifndef UTILITIES_DECLARATIONS
#include "utilities.hpp"
#endif


namespace surrogate {


FLOAT determinationCoeff (
  const std::vector<FLOAT>& responses, 
  const std::vector<FLOAT>& approximations
);


} // namespace surrogate


#endif // REGRESSION_DECLARATIONS

