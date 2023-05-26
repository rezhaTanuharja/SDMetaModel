/**
  * @file MassSpringDamper_DMC.cpp 
  *
  * @brief
  * implements Direct MC method for stochastic mass-spring-damper model 
  *
  * @author
  * Tanuharja R.A. @n 
  * contact: rezha.tanuharja@tum.de 
  */

#include "MassSpringDamper_DMC.hpp" 

namespace Stochastic::DirectMC {


MassSpringDamper::MassSpringDamper ( const Z nDOFs, const Z Dim ) :
  nDOFs_ ( nDOFs ), Dim_ ( Dim ), DampingRatio_ ( 0.0 )
{}


void MassSpringDamper::SetMassCoeffs ( const VectorR& MassCoeffs ) {


  #ifdef DEBUG_RSMSD 

  if ( MassCoeffs.size() != nDOFs_ * Dim_ ) {
    throw std::runtime_error (
      "MSD_DMC::SetMassCoeffs: incorrect size"
    );
  }
  #endif // DEBUG_RSMSD 
  
  
  MassCoeffs_ = MassCoeffs; 
}

void MassSpringDamper::SetSpringCoeffs ( const VectorR& SpringCoeffs ) {


  #ifdef DEBUG_RSMSD 

  if ( SpringCoeffs.size() != nDOFs_ * Dim_ ) {
    throw std::runtime_error (
      "MSD_DMC::SetSpringCoeffs: incorrect size"
    );
  }
  #endif // DEBUG_RSMSD 
  
  
  SpringCoeffs_ = SpringCoeffs; 
}


void MassSpringDamper::SetForceCoeffs ( const VectorC& ForceCoeffs ) {


  #ifdef DEBUG_RSMSD 

  if ( ForceCoeffs.size() != nDOFs_ * Dim_ ) {
    throw std::runtime_error (
      "MSD_DMC::SetForceCoeffs: incorrect size"
    );
  }
  #endif // DEBUG_RSMSD 
  
  
  ForceCoeffs_ = ForceCoeffs; 
}


void MassSpringDamper::SetDampingRatio ( const R DampingRatio ) {


  #ifdef DEBUG_RSMSD 

  if ( DampingRatio < 0 ) {
    throw std::runtime_error (
      "MSD_DMC::SetDampingRatio: negative DampingRatio"
    );
  }

  #endif 
  
  
  DampingRatio_ = DampingRatio; 
}


VectorC MassSpringDamper::ComputeResponses (
  const VectorR& RandomBasis, 
  const R Omega 
) {

  auto nPoints = RandomBasis.size() / Dim_; 


  #ifdef DEBUG_RSMSD 

  if ( RandomBasis.size() != nPoints * Dim_ ){
    throw std::runtime_error (
      "MSD_DMC::ComputeResponses: RandomBasis size not a multiple of Dim_"
    );
  }

  #endif // DEBUG_RSMSD 
  
  
  Eigen::Map<const MatrixXR> randomBasis (
    RandomBasis.data(), Dim_, nPoints 
  );


  // Compute realisations of masses -----------------------------------

  VectorR RandomMasses ( nDOFs_ * nPoints ); 

  Eigen::Map<MatrixXR> randomMasses (
    RandomMasses.data(), nDOFs_, nPoints 
  );

  Eigen::Map<const MatrixXR> massCoeffs (
    MassCoeffs_.data(), nDOFs_, Dim_ 
  );

  randomMasses = massCoeffs * randomBasis; 


  // Compute realisations of springs ----------------------------------

  VectorR RandomSprings ( nDOFs_ * nPoints ); 

  Eigen::Map<MatrixXR> randomSprings (
    RandomSprings.data(), nDOFs_, nPoints 
  );

  Eigen::Map<const MatrixXR> springCoeffs (
    SpringCoeffs_.data(), nDOFs_, Dim_ 
  );

  randomSprings = springCoeffs * randomBasis; 


  // Compute realisations of forces -----------------------------------
  
  VectorC RandomForces ( nDOFs_ * nPoints ); 

  Eigen::Map<MatrixXC> randomForces (
    RandomForces.data(), nDOFs_, nPoints 
  );

  Eigen::Map<const MatrixXC> forceCoeffs (
    ForceCoeffs_.data(), nDOFs_, Dim_ 
  );

  randomForces = forceCoeffs * randomBasis; 


  // Compute realisations of responses --------------------------------
  
  VectorC Result ( nDOFs_ * nPoints ); 

  Eigen::Map<MatrixXC> result (
    Result.data(), nDOFs_, nPoints 
  );

  namespace MSD = SimpleModels::MassSpringDamper; 

  for ( auto i = 0; i < nPoints; i++ ) {

    VectorR Masses ( 
      RandomMasses.begin() + i * nDOFs_, 
      RandomMasses.begin() + i * nDOFs_ + nDOFs_ 
    ); 

    VectorR Springs ( 
      RandomSprings.begin() + i * nDOFs_, 
      RandomSprings.begin() + i * nDOFs_ + nDOFs_ 
    ); 

    VectorC Forces (
      RandomForces.begin() + i * nDOFs_, 
      RandomForces.begin() + i * nDOFs_ + nDOFs_ 
    );

    auto MassMatrix = MSD::MassMatrix<R> (
      Masses
    );

    auto StiffnessMatrix = MSD::StiffnessMatrix<R> (
      Springs
    );

    auto DynamicStiffnessMatrix = MSD::DynamicStiffnessMatrix<R,C> (
      MassMatrix, StiffnessMatrix, DampingRatio_, Omega 
    );

    Eigen::Map<const MatrixXC> dynamicStiffnessMatrix (
      DynamicStiffnessMatrix.data(), nDOFs_, nDOFs_ 
    );

    Eigen::Map<const VectorXC> forces (
      Forces.data(), nDOFs_ 
    );

    result.col(i) = dynamicStiffnessMatrix.partialPivLu().solve(forces);

  }

  return Result; 

} // MassSpringDamper::ComputeResponses 


} // namespace Stochastic::DirectMC 

