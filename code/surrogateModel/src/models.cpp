/**
 * @file models.cpp
 *
 * @brief
 * Contains implementations of surrogate models.
 *
 * @author Rezha Adrian Tanuharja <rezhadr@outlook.com>
 * @date 2023-11-02
*/

#ifndef MODELS_DECLARATIONS
#include "models.hpp"
#endif


namespace surrogate {


RPCE::RPCE ( const INT nDOFs, const INT dimension ) :
  nDOFs_ ( nDOFs ), dimension_ ( dimension ) 
{}


void RPCE::setNumIndices ( const std::vector<INT>& indices ) {
  numIndices_ = indices;
}


void RPCE::setDenIndices ( const std::vector<INT>& indices ) {
  denIndices_ = indices;
}


void RPCE::train (
  const std::vector<FLOAT>& args, 
  const std::vector<FLOAT>& responses 
) {

  INT nPoints = args.size() / dimension_;

  INT nNumBasis = numIndices_.size() / dimension_;
  INT nDenBasis = denIndices_.size() / dimension_;

  std::vector<FLOAT> numBasis = hermiteBasis ( numIndices_, args, dimension_ );
  std::vector<FLOAT> denBasis = hermiteBasis ( denIndices_, args, dimension_ );

  Eigen::Map<MatrixXF> numPsi (
    numBasis.data(), nNumBasis, nPoints
  );
  Eigen::Map<MatrixXF> denPsi (
    denBasis.data(), nDenBasis, nPoints
  );

  Eigen::Map<const MatrixXF> result (
    responses.data(), nDOFs_, nPoints
  );

  numCoeffs_ = std::vector<FLOAT> ( nDOFs_ * nNumBasis, 0.0 );
  denCoeffs_ = std::vector<FLOAT> ( nDOFs_ * nDenBasis, 0.0 );

  #pragma omp parallel for
  for ( INT i = 0; i < nDOFs_; i++ ) {

    DiagMatrixXF M  = result.row(i).asDiagonal();
    DiagMatrixXF Mh = result.row(i).conjugate().asDiagonal();

    MatrixXF Aupper ( nNumBasis, nNumBasis + nDenBasis );

    Aupper << 
      (  numPsi   *   numPsi.transpose() ),
      ( -numPsi * M * denPsi.transpose() );

    MatrixXF Alower ( nDenBasis, nNumBasis + nDenBasis );

    Alower << 
      ( -denPsi * Mh   *   numPsi.transpose() ),
      (  denPsi * Mh * M * denPsi.transpose() );

    MatrixXF A ( nNumBasis + nDenBasis, nNumBasis + nDenBasis );
    A << Aupper, Alower;

    MatrixXF V = A.bdcSvd ( Eigen::ComputeFullV ).matrixV();

    VectorXF r = V.col( V.cols() - 1 );

    for ( INT j = 0; j < nNumBasis; j++ ) {
      numCoeffs_[ i * nNumBasis + j ] = r(j);
    }

    for ( INT j = 0; j < nDenBasis; j++ ) {
      denCoeffs_[ i * nDenBasis + j ] = r(nNumBasis + j);
    }

  }

} // RPCE::train 


void RPCE::trainSparse (
  const std::vector<FLOAT>& args,
  const std::vector<FLOAT>& responses, 
  const FLOAT epsilon
) {



} // RPCE::trainSparse


std::vector<FLOAT> RPCE::computeResponse ( 
  const std::vector<FLOAT>& args 
) const {

  INT nPoints = args.size() / dimension_;

  INT nNumBasis = numIndices_.size() / dimension_;
  INT nDenBasis = denIndices_.size() / dimension_;

  std::vector<FLOAT> numBasis = hermiteBasis ( numIndices_, args, dimension_ );
  std::vector<FLOAT> denBasis = hermiteBasis ( denIndices_, args, dimension_ );

  Eigen::Map<MatrixXF> numPsi ( numBasis.data(), nNumBasis, nPoints );
  Eigen::Map<MatrixXF> denPsi ( denBasis.data(), nDenBasis, nPoints );

  Eigen::Map<const MatrixXF> numCoeffs ( numCoeffs_.data(), nNumBasis, nDOFs_ );
  Eigen::Map<const MatrixXF> denCoeffs ( denCoeffs_.data(), nDenBasis, nDOFs_ );

  std::vector<FLOAT> Responses ( nDOFs_ * nPoints );

  Eigen::Map<MatrixXF> responses ( Responses.data(), nDOFs_, nPoints );

  responses = 
    ( numCoeffs.transpose() * numPsi ).array() / 
    ( denCoeffs.transpose() * denPsi ).array();

  return Responses;

} // surrogate::RPCE::computeResponse 


} // namespace surrogate

