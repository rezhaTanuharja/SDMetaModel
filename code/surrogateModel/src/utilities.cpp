/**
 * @file utilities.cpp
 *
 * @brief
 * Contains implementations of various useful functions.
 *
 * @author Rezha Adrian Tanuharja <rezhadr@outlook.com>
 * @date 2023-11-02
*/

#ifndef UTILITIES_DECLARATIONS
#include "utilities.hpp"
#endif


namespace surrogate {


INT factorial ( const INT n ) {

  if ( n == 0 ) {
    return 1;
  }

  return n * factorial ( n - 1 );

} // factorial


INT binomialCoeff ( const INT n, const INT k ) {

  if ( n < k ) { return 0; }

  if ( k == 0 || k == n ) { return 1; }

  return 
    binomialCoeff ( n - 1, k - 1 ) +
    binomialCoeff ( n - 1, k     );

} // binomialCoeff


FLOAT hermiteFunction ( const INT index, const FLOAT x ) {

  if ( index == 0 ) { return 1.0; }
  if ( index == 1 ) { return x; }

  return (
    x             * hermiteFunction ( index - 1, x ) - 
    ( index - 1 ) * hermiteFunction ( index - 2, x )
  );

} // hermiteFunction


FLOAT normalVariable ( const FLOAT mean, const FLOAT stdDev, const INT seed ) {

  std::mt19937 generator ( seed );

  std::normal_distribution<FLOAT> normalDist ( mean, stdDev );

  return normalDist ( generator );

} // normalVariable


std::vector<INT> randomSubset (
  const std::vector<INT>& set, const INT size, const INT seed
) {

  std::mt19937 generator ( seed );

  std::vector<INT> order ( set.size() );
  std::copy (
    set.begin(), set.end(), order.begin()
  );

  std::shuffle ( order.begin(), order.end(), generator );

  std::vector<INT> subset;

  for ( INT num : order ) {
    if ( subset.size() == size ) {
      break;
    }
    subset.push_back ( num );
  }

  std::sort ( subset.begin(), subset.end() );

  return subset;

} // randomSubset


void saveVector (
  const std::string& filename, const std::vector<FLOAT>& vector
) {

  std::ofstream file ( filename, std::ios::binary | std::ios::app );

  if ( file.is_open() ) {

    INT size = vector.size();
    file.write ( reinterpret_cast<const char*> ( &size ), sizeof ( size ) );

    file.write (
      reinterpret_cast<const char*> ( vector.data() ),
      vector.size() * sizeof ( FLOAT )
    );

    file.close();

  } // file.is_open()

} // saveVector


void loadVector (
  const std::string& filename, std::vector<FLOAT>& vector
) {

  std::ifstream file ( filename, std::ios::binary );

  if ( file.is_open() ) {

    vector.clear();

    while ( true ) {

      INT size; 
      if ( !file.read( reinterpret_cast<char*>(&size), sizeof(size) )) {
        break;
      }

      std::vector<FLOAT> data ( size );
      if ( 
        !file.read( reinterpret_cast<char*>(data.data()), size * sizeof(FLOAT) )
      ) {
        break;
      }

      vector.insert ( vector.end(), data.begin(), data.end() );

    }

    file.close();

  } // file.is_open()

} // loadVector


void loadVector (
  const std::string& filename, std::vector<INT>& vector
) {

    std::ifstream inputFile(filename);

    if (!inputFile.is_open()) {
        return;
    }

    INT value;
    while (inputFile >> value) {
        vector.push_back(value);
    }

    // Close the file
    inputFile.close();


} // loadVector


} // namespace surrogate

