/**
  * @file NIRPCE.hpp 
  *
  * @brief
  * declares supports for NIRPCE modelling 
  *
  * @author
  * Tanuharja R.A. @n 
  * contact: rezha.tanuharja@tum.de 
  */

#ifndef NIRPCE_DECLARATIONS 
#define NIRPCE_DECLARATIONS 

#ifndef STOCHASTIC_DECLARATIONS 
  #include "Stochastic.hpp" 
#endif 

namespace Stochastic::Surrogate {


/**
  * @class NIRPCE 
  *
  * @brief 
  * Model response as ratio of two PCEs 
  * Implemented in @ref _NIRPCE_cpp_ 
  */
class NIRPCE {

  VectorZ NumIndices_; 
  VectorZ DenIndices_; 

  VectorC NumCoeffs_; 
  VectorC DenCoeffs_; 

  [[maybe_unused]] Z nDOFs_; 
  [[maybe_unused]] Z Dim_;


  public: 


  /**
    * @brief 
    * Constructs NIRPCE model for mass-spring-damper model. @n 
    * Implemented in @ref _NIRPCE_cpp_  
    *
    * @param nDOFs number of degree of freedoms 
    * @param Dim   number of independent random variables 
    */
  NIRPCE ( const Z nDOFs, const Z Dim ); 


  /**
    * @brief 
    * Set indices of basis functions for numerator 
    * Implemented in @ref _NIRPCE_cpp_  
    *
    * @param NumIndices indices of basis functions 
    */
  void SetNumIndices ( const VectorZ& NumIndices ); 


  /**
    * @brief 
    * Set indices of basis functions for denominator  
    * Implemented in @ref _NIRPCE_cpp_  
    *
    * @param NumIndices indices of basis functions 
    */
  void SetDenIndices ( const VectorZ& DenIndices ); 


  /**
    * @brief 
    * Compute numerator and denominator coefficients 
    * Implemented in @ref _NIRPCE_cpp_  
    *
    * @param InputVars realisations of input random variables 
    * @param Responses outputs associated with input random variables 
    */
  void Train (
    const VectorC& InputVars, 
    const VectorC& Responses 
  );


  /**
    * @brief 
    * Compute surrogate model response for given random input variables
    * Implemented in @ref _NIRPCE_cpp_  
    *
    * @param InputVars realisations of input random variables 
    *
    * @return (predicted) outputs associated with input random variables 
    */
  VectorC ComputeResponse ( const VectorC& InputVars ); 


  #ifdef DEBUG_RSMSD 

  VectorC NumCoeffs () const {
    return NumCoeffs_; 
  }

  VectorC DenCoeffs () const {
    return DenCoeffs_; 
  }

  #endif // DEBUG_RSMSD 


}; // class NIRPCE 


} // namespace Stochastic::Surrogate 

#endif // NIRPCE_DECLARATIONS 

