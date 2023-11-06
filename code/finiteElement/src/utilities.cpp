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


namespace mFEM {


IntegrationPoints gaussQuadrature ( const INT n ) {

  #ifdef ADDITIONAL_DEBUG

  if ( n > 3 || n < 0 ) {
    throw std::runtime_error (
      "mFEM::gaussQuadrature: do not support n > 3 or n < 0!"
    );
  }

  #endif // ADDITIONAL_DEBUG

  std::vector<FLOAT> points;
  std::vector<FLOAT> weights;

  if ( n == 3 ) {
    points  = { -std::sqrt(0.6), 0.0, std::sqrt(0.6) };
    weights = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

    return { points, weights };
  }

  if ( n == 2 ) {
    points  = { -1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0) };
    weights = { 1.0, 1.0 };

    return { points, weights };
  }

  points  = { 0.0 };
  weights = { 2.0 };

  return { points, weights };

} // gaussQuadrature


FLOAT normalVariable ( const FLOAT mean, const FLOAT stdDev, const INT seed ) {

  std::mt19937 generator ( seed );

  std::normal_distribution<FLOAT> normalDist ( mean, stdDev );

  return normalDist ( generator );

} // normalVariable


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


void saveVector (
  const std::string& filename, const std::vector<INT>& vector
) {

  std::ofstream file ( filename, std::ios::app );

  if ( file.is_open() ) {

    for ( const INT value : vector ) {
      file << value << std::endl;
    }

    file.close();

  } // file.is_open()

} // saveVector


void loadVector (
  const std::string& filename, std::vector<INT>& vector
) {

  std::ifstream file ( filename );

  if ( file.is_open() ) {

    INT value;
    while ( file >> value ) {
      vector.push_back ( value );
    }

    file.close();

  } // file.is_open()

} // loadVector


MatrixXF zeroPadRows (
  const MatrixXF& matrix,
  const std::vector<INT>& internalDofIds, 
  const std::vector<INT>& boundaryDofIds
) {

  INT nInt = internalDofIds.size();
  INT nExt = boundaryDofIds.size();

  INT nRows = internalDofIds.size() + boundaryDofIds.size();

  INT p = 0;
  INT q = 0;

  MatrixXF paddedMatrix ( nRows, matrix.cols() );

  for ( INT i = 0; i < nRows; i++ ) {

    if ( p >= nInt ) {
      break;
    }

    if ( q >= nExt ) {
      paddedMatrix.row(i) = matrix.row(p);
      p++;
      continue;
    }

    if ( internalDofIds[p] < boundaryDofIds[q] ) {

      paddedMatrix.row(i) = matrix.row(p);
      p++;
      continue;

    }

    q++;

  }


  return paddedMatrix;

} // zeroPadRows


VectorXF combine (
  VectorXF Xn, VectorXF Xb, 
  const std::vector<INT>& internalDofIds, 
  const std::vector<INT>& boundaryDofIds
) {

  VectorXF X ( internalDofIds.size() + boundaryDofIds.size() );

  for ( INT i = 0; i < internalDofIds.size(); i++ ) {
    X.coeffRef ( internalDofIds[i] ) = Xn[i];
  }

  for ( INT i = 0; i < boundaryDofIds.size(); i++ ) {
    X.coeffRef ( boundaryDofIds[i] ) = Xb[i];
  }

  return X;

} // combine 


} // namespace mFEM

