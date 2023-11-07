/**
 * @file basisFunctions.cpp
 *
 * @brief
 * Contains implementations of basis functions for PCE and RPCE.
 *
 * @author Rezha Adrian Tanuharja <rezhadr@outlook.com>
 * @date 2023-11-02
*/

#ifndef BASIS_FUNCTIONS_DECLARATIONS
#include "basisFunctions.hpp"
#endif


namespace surrogate {


/**
 * @private
 * A subroutine to help generate indices.
*/
void generateNextIndices (
  std::vector<INT>& indices, std::vector<bool>& lastSet, const INT n 
) {

  std::vector<INT> nextIndices;

  std::vector<INT> integers ( n );
  std::iota ( integers.begin(), integers.end(), 1 );

  std::next_permutation ( lastSet.begin(), lastSet.end() );

  INT m = 0; 
  for ( const bool& active : lastSet ) {
    if ( active ) {
      nextIndices.push_back ( integers[m] );
    }
    m++;
  }

  nextIndices.push_back ( n + 1 );

  std::adjacent_difference (
    nextIndices.begin(), nextIndices.end(), 
    nextIndices.begin()
  );

  std::transform (
    nextIndices.begin(), nextIndices.end(), 
    nextIndices.begin(),
    [=]( const INT m ) { return m - 1; }
  );

  indices.insert ( indices.end(), nextIndices.begin(), nextIndices.end() );

} // generateNextIndices


std::vector<INT> generateIndices (
  const INT dimension, const INT maxSum
) {

  std::vector<INT> indices;
  indices.reserve (
    dimension * binomialCoeff ( maxSum + dimension, dimension )
  );

  for ( INT i = 0; i < dimension; i++ ) {
    indices.push_back ( 0 );
  }

  for ( INT i = 0; i < maxSum; i++ ) {
    
    INT jMax = binomialCoeff ( dimension + i, dimension - 1 );

    std::vector<bool> lastSet ( dimension + i );

    std::transform (
      lastSet.begin(), lastSet.begin() + dimension - 1,
      lastSet.begin(),
      [=]( const auto m ) { return true; }
    );

    for ( INT j = 0; j < jMax; j++ ) {
      generateNextIndices ( indices, lastSet, dimension + i );
    }

  }

  return indices;

} // generateIndices


std::vector<FLOAT> hermiteBasis (
  const std::vector<INT>& indices, 
  const std::vector<FLOAT>& args, 
  const INT dimension
) {

  INT nProducts = indices.size() / dimension;
  INT nSamples  = args.size() / dimension;

  std::vector<FLOAT> basis;
  basis.reserve ( nProducts * nSamples );

  for ( INT i = 0; i < nSamples; i++ ) {
    for ( INT j = 0; j < nProducts; j++ ) {

      FLOAT unity = 1.0;

      basis.push_back (
        std::transform_reduce (
          indices.begin() + j * dimension, 
          indices.begin() + j * dimension + dimension, 
          args.begin() + i * dimension, 

          unity, 

          [](const FLOAT N, const FLOAT P ) { return N * P; }, 

          [](const INT idx, const FLOAT x ) {
            return hermiteFunction ( idx, x ) /
            std::sqrt ( factorial ( idx ) );
          }

        )
      );

    }
  }

  return basis;

} // hermiteBasis


} // namespace surrogate

