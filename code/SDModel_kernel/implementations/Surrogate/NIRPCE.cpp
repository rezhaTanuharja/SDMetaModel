/**
  * @file NIRPCE.cpp 
  *
  * @brief
  * implements supports for NIRPCE modelling 
  *
  * @author
  * Tanuharja R.A. @n 
  * contact: rezha.tanuharja@tum.de 
  */

#include "NIRPCE.hpp" 

namespace Stochastic::Surrogate {


NIRPCE::NIRPCE ( const Z nDOFs, const Z Dim ) :
  nDOFs_ ( nDOFs ), Dim_ ( Dim ) 
{}


void NIRPCE::SetNumIndices ( const VectorZ& NumIndices ) {
  NumIndices_ = NumIndices; 
}


void NIRPCE::SetDenIndices ( const VectorZ& DenIndices ) {
  DenIndices_ = DenIndices; 
}


void NIRPCE::Train ( const VectorC& InputVars, const VectorC& Responses ) {

  auto nPoints = InputVars.size() / Dim_; 

  auto nNumBasis = NumIndices_.size() / Dim_; 
  auto nDenBasis = DenIndices_.size() / Dim_; 

  auto NumBasis = BasisFunctions::ComputeHermiteBasis<Z,R,C> (
    NumIndices_, InputVars, Dim_ 
  );

  auto DenBasis = BasisFunctions::ComputeHermiteBasis<Z,R,C> (
    DenIndices_, InputVars, Dim_ 
  );

  VectorC Coeffs ( nDOFs_ * ( nNumBasis + nDenBasis ) );

  Eigen::Map<MatrixXC> NumPsi (
    NumBasis.data(), nNumBasis, nPoints 
  );

  Eigen::Map<MatrixXC> DenPsi (
    DenBasis.data(), nDenBasis, nPoints 
  );

  Eigen::Map<MatrixXC> coeffs (
    Coeffs.data(), nNumBasis + nDenBasis, nDOFs_ 
  );

  Eigen::Map<const MatrixXC> Result (
    Responses.data(), nDOFs_, nPoints 
  );

  NumCoeffs_ = VectorC ( nDOFs_ * nNumBasis, 0.0 );
  DenCoeffs_ = VectorC ( nDOFs_ * nDenBasis, 0.0 );

  for ( auto i = 0; i < nDOFs_; i++ ) {

    DiagMatrixXC M  = Result.row(i).asDiagonal(); 
    DiagMatrixXC Mh = Result.row(i).conjugate().asDiagonal(); 

    MatrixXC Aupper ( nNumBasis, nNumBasis + nDenBasis );
    Aupper << (  NumPsi   *   NumPsi.transpose() ), 
              ( -NumPsi * M * DenPsi.transpose() );

    MatrixXC Alower ( nDenBasis, nNumBasis + nDenBasis ); 
    Alower << ( -DenPsi * Mh   *   NumPsi.transpose() ),
              (  DenPsi * Mh * M * DenPsi.transpose() );

    MatrixXC A ( nNumBasis + nDenBasis, nNumBasis + nDenBasis ); 
    A << Aupper, Alower; 

    MatrixXC V = A.bdcSvd ( Eigen::ComputeFullV ).matrixV();

    VectorXC r = V.col( V.cols() - 1 );

    for ( auto j = 0; j < nNumBasis; j++ ) {
      NumCoeffs_[ i * nNumBasis + j ] = r(j);
    }

    for ( auto j = 0; j < nDenBasis; j++ ) {
      DenCoeffs_[ i * nDenBasis + j ] = r( nNumBasis + j );
    }

  } // compute coeffs for each DOF 

} // NIRPCE::Train 


VectorC NIRPCE::ComputeResponse ( const VectorC& InputVars ) {

  auto nPoints = InputVars.size() / Dim_; 

  auto nNumBasis = NumIndices_.size() / Dim_;
  auto nDenBasis = DenIndices_.size() / Dim_;

  auto NumBasis = BasisFunctions::ComputeHermiteBasis<Z,R,C> (
    NumIndices_, InputVars, Dim_ 
  );

  auto DenBasis = BasisFunctions::ComputeHermiteBasis<Z,R,C> (
    DenIndices_, InputVars, Dim_ 
  );

  Eigen::Map<MatrixXC> NumPsi (
    NumBasis.data(), nNumBasis, nPoints 
  );

  Eigen::Map<MatrixXC> DenPsi (
    DenBasis.data(), nDenBasis, nPoints 
  );

  Eigen::Map<MatrixXC> numCoeffs (
    NumCoeffs_.data(), nNumBasis, nDOFs_ 
  );

  Eigen::Map<MatrixXC> denCoeffs (
    DenCoeffs_.data(), nDenBasis, nDOFs_ 
  );

  VectorC Result ( nDOFs_ * nPoints ); 

  Eigen::Map<MatrixXC> result (
    Result.data(), nDOFs_, nPoints 
  );

  result = ( numCoeffs.transpose() * NumPsi ).array() /
           ( denCoeffs.transpose() * DenPsi ).array();

  return Result; 

}



} // namespace Stochastic::Surrogate 

