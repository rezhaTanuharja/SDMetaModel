/**
  * @file MassSpringDamper_DMC.hpp 
  *
  * @brief
  * declares Direct MC method for stochastic mass-spring-damper model 
  *
  * @author
  * Tanuharja R.A. @n 
  * contact: rezha.tanuharja@tum.de 
  */

#ifndef MASS_SPRING_DAMPER_DMC_DECLARATIONS 
#define MASS_SPRING_DAMPER_DMC_DECLARATIONS 

#ifndef STOCHASTIC_DECLARATIONS 
  #include "Stochastic.hpp" 
#endif 

namespace Stochastic::DirectMC {


/**
  * @class MassSpringDamper
  *
  * @brief 
  * Supports direct MC simulation of stochastic mass-spring-damper model. @n 
  * Implemented in @ref _MassSpringDamper_DMC_cpp_ 
  */
class MassSpringDamper {

  VectorR MassCoeffs_; 
  VectorR SpringCoeffs_; 

  VectorC ForceCoeffs_; 

  Z nDOFs_; 
  Z Dim_; 

  R DampingRatio_; 


  public: 


  /**
    * @brief 
    * Constructs Direct MC object for mass-spring-damper model. 
    * Constructs undamped model by default. @n 
    * Implemented in @ref _MassSpringDamper_DMC_cpp_ 
    *
    * @param nDOFs number of degree of freedoms 
    * @param Dim   number of independent random variables 
    */
  MassSpringDamper ( const Z nDOFs, const Z Dim );


  /**
    * @brief 
    * Set mass coeffs of stochastic basis functions. @n
    * Implemented in @ref _MassSpringDamper_DMC_cpp_ 
    *
    * @param MassCoeffs nDOFs x Dim matrix in column major storage order 
    */
  void SetMassCoeffs ( const VectorR& MassCoeffs ); 


  /**
    * @brief 
    * Set spring coeffs of stochastic basis functions. @n
    * Implemented in @ref _MassSpringDamper_DMC_cpp_ 
    *
    * @param SpringCoeffs nDOFs x Dim matrix in column major storage order 
    */
  void SetSpringCoeffs ( const VectorR& SpringCoeffs ); 


  /**
    * @brief 
    * Set force coeffs of stochastic basis functions. @n
    * Implemented in @ref _MassSpringDamper_DMC_cpp_ 
    *
    * @param ForceCoeffs nDOFs x Dim matrix in column major storage order 
    */
  void SetForceCoeffs ( const VectorC& ForceCoeffs ); 


  /**
    * @brief 
    * Set DampingRatio for damped model. @n 
    * Implemented in @ref _MassSpringDamper_DMC_cpp_ 
    *
    * @param DampingRatio structural damping parameter 
    */
  void SetDampingRatio ( const R DampingRatio ); 


  /**
    * @brief 
    * Compute responses for given realisations of random variables. @n
    * Implemented in @ref _MassSpringDamper_DMC_cpp_ 
    *
    * @param RandomBasis realisation of random basis functions 
    * @param Omega       angular speed 
    *
    * @return nDOFs x nPoints matrix in column major storage order 
    */
  VectorC ComputeResponses (
    const VectorR& RandomBasis, 
    const R Omega 
  );


  #ifdef DEBUG_RSMSD 

  VectorR MassCoeffs () const {
    return MassCoeffs_; 
  }

  VectorR SpringCoeffs () const {
    return SpringCoeffs_; 
  }
 
  VectorC ForceCoeffs () const {
    return ForceCoeffs_; 
  }

  R DampingRatio () const {
    return DampingRatio_; 
  }
 
  #endif 
  
  
}; // class MassSpringDamper 


} // namespace Stochastic::DirectMC 


#endif // MASS_SPRING_DAMPER_DMC_DECLARATIONS 

