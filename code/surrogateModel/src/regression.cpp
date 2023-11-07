/**
 * @file regression.cpp
 *
 * @brief
 * Contains implementations of various functions for regression.
 *
 * @author Rezha Adrian Tanuharja <rezhadr@outlook.com>
 * @date 2023-11-02
*/

#ifndef REGRESSION_DECLARATIONS
#include "regression.hpp"
#endif


namespace surrogate {


FLOAT determinationCoeff (
  const std::vector<FLOAT>& responses, 
  const std::vector<FLOAT>& approximations
) {

  INT n = responses.size();

  FLOAT zero = 0.0;

  FLOAT I = std::transform_reduce (
    responses.begin(), responses.end(), 
    approximations.begin(), 

    zero, 

    [](const FLOAT N, const FLOAT P) { return N + P; },

    [](const FLOAT p, const FLOAT q) {
      return ( p - q ) * ( p - q );
    }

  );

  I /= n;

  FLOAT mean = std::accumulate ( responses.begin(), responses.end(), 0.0 );
  mean /= n;

  std::vector<FLOAT> offsets ( n );

  std::transform (
    responses.begin(), responses.end(), 
    offsets.begin(), 

    [=](const auto m) {
      return ( m - mean ) * ( m - mean );
    }
  );

  FLOAT var = std::accumulate ( offsets.begin(), offsets.end(), 0.0 );
  var /= n - 1;

  return 1 - I / var;

} // determinationCoeff


} // namespace surrogate

