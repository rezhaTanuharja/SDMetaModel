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

#ifndef REGRESSION_DECLARATIONS
#include "regression.hpp"
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

  INT nPoints = args.size() / dimension_;

  // Compute reference coefficient of determination i.e., a baseline

  std::vector<INT> pointID ( nPoints );
  std::iota( pointID.begin(), pointID.end(), 0);

  std::vector<INT> refPoints = surrogate::randomSubset( 
    pointID, nPoints, 0
  );

  std::vector<FLOAT> refArgs;
  std::vector<FLOAT> refResp;

  for ( INT i = 0; i < nPoints; i++ ) {

    if ( i >= INT ( nPoints * 0.95 ) ) {
      break;
    }

    for ( INT j = 0; j < dimension_; j++ ) {
      refArgs.push_back(
        args[
          dimension_ * refPoints[i] + j
        ]
      );
    }

    for ( INT j = 0; j < nDOFs_; j++ ) {
      refResp.push_back(
        responses[
          nDOFs_ * refPoints[i] + j
        ]
      );
    }

  }

  this -> train( refArgs, refResp );
  FLOAT eps = determinationCoeff ( responses, this -> computeResponse ( args ) );

  // Backward step for numerator

  INT nNumBasis = numIndices_.size() / dimension_;

  for ( INT i = 0; i < nNumBasis; i++ ) {

    std::vector<INT> sortedPoint = surrogate::randomSubset( 
      pointID, nPoints, i
    );

    std::vector<FLOAT> trainArgs;
    std::vector<FLOAT> trainResp;

    for ( INT j = 0; j < nPoints; j++ ) {

      if ( j >= INT ( nPoints * 0.95 ) ) {
        break;
      }

      for ( INT k = 0; k < dimension_; k++ ) {
        trainArgs.push_back(
          args[
            dimension_ * sortedPoint[j] + k
          ]
        );
      }

      for ( INT k = 0; k < nDOFs_; k++ ) {
        trainResp.push_back(
          responses[
            nDOFs_ * sortedPoint[j] + k
          ]
        );
      }

    }

    INT n = nNumBasis - 1 - i;

    std::vector<INT> oldIndices ( numIndices_.size() );
    std::copy ( numIndices_.begin(), numIndices_.end(), oldIndices.begin() );

    std::vector<INT> newIndices ( numIndices_.size() );
    std::copy ( numIndices_.begin(), numIndices_.end(), newIndices.begin() );

    newIndices.erase( 
      newIndices.begin() + n * dimension_,
      newIndices.begin() + n * dimension_ + dimension_
    );

    this -> setNumIndices( newIndices ) ;
    this -> train( trainArgs, trainResp );

    FLOAT newEps = determinationCoeff ( 
      responses, this -> computeResponse ( args) 
    );

    if ( newEps > eps ) {
      eps = newEps;
    }
    else if ( eps - newEps > epsilon * ( 1.0 - eps ) ) {
      this -> setNumIndices( oldIndices );
    }
    else {
      eps = newEps;
    }

  }

  // Backward step for denominator

  INT nDenBasis = denIndices_.size() / dimension_;

  for ( INT i = 0; i < nDenBasis; i++ ) {

    std::vector<INT> sortedPoint = surrogate::randomSubset( 
      pointID, nPoints,i + 5
    );

    std::vector<FLOAT> trainArgs;
    std::vector<FLOAT> trainResp;

    for ( INT j = 0; j < nPoints; j++ ) {

      if ( j >= INT ( nPoints * 0.95 ) ) {
        break;
      }

      for ( INT k = 0; k < dimension_; k++ ) {
        trainArgs.push_back(
          args[
            dimension_ * sortedPoint[j] + k
          ]
        );
      }

      for ( INT k = 0; k < nDOFs_; k++ ) {
        trainResp.push_back(
          responses[
            nDOFs_ * sortedPoint[j] + k
          ]
        );
      }

    }

    INT n = nDenBasis - 1 - i;

    std::vector<INT> oldIndices ( denIndices_.size() );
    std::copy ( denIndices_.begin(), denIndices_.end(), oldIndices.begin() );

    std::vector<INT> newIndices ( denIndices_.size() );
    std::copy ( denIndices_.begin(), denIndices_.end(), newIndices.begin() );

    newIndices.erase( 
      newIndices.begin() + n * dimension_,
      newIndices.begin() + n * dimension_ + dimension_
    );

    this -> setDenIndices( newIndices ) ;
    this -> train( trainArgs, trainResp );

    FLOAT newEps = determinationCoeff ( 
      responses, this -> computeResponse ( args ) 
    );

    if ( newEps > eps ) {
      eps = newEps;
    }
    else if ( eps - newEps > epsilon * ( 1.0 - eps ) ) {
      this -> setDenIndices( oldIndices );
    }
    else {
      eps = newEps;
    }

  }

  this -> train( args, responses );

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

