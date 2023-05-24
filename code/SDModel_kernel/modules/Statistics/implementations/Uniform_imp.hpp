/**
  * @file Uniform_imp.hpp 
  *
  * @brief
  * implements template functions to generate RVs from uniform distributions 
  *
  * @anchor _Uniform_imp_hpp_ 
  *
  * @author
  * Tanuharja R.A. @n 
  * contact: rezha.tanuharja@tum.de 
  */

#ifndef UNIFORM_IMPLEMENTATIONS 
#define UNIFORM_IMPLEMENTATIONS 

#ifndef UNIFORM_DECLARATIONS 
  #include "Uniform.hpp" 
#endif 

namespace Statistics {

namespace Uniform {


template < typename Z, typename R, typename C >
Vector<C> RandomSampling (
  const Z nPoints, const Z Dim, 
  const R LowerBound, const R UpperBound
) {
  

  #ifdef DEBUG_STATISTICS 

  if ( nPoints < 1 || nPoints > 80000 ) {
    throw std::runtime_error (
      "Uniform::RandomSampling: nPoints not in [1,80000] range"
    );
  }

  if ( Dim < 1 || Dim > 100 ) {
    throw std::runtime_error (
      "Uniform::RandomSampling: Dim not in [1,100] range"
    );
  }

  if ( !( LowerBound < UpperBound ) ) {
    throw std::runtime_error (
      "Uniform::RandomSampling: LowerBound >= UpperBound"
    );
  }

  #endif // DEBUG_STATISTICS 
  
  
  std::random_device device; 
  std::mt19937 generator ( device () );

  std::uniform_real_distribution<R> RandomVariable ( LowerBound, UpperBound );

  Vector<C> result ( nPoints * Dim ); 

  for ( auto i = 0; i < result.size(); i++ ) {
    result[i] = C ( RandomVariable ( generator ) );
  }

  return result;

} // RandomSampling 


} // namespace Uniform 

} // namespace Statistics 

#endif // UNIFORM_IMPLEMENTATIONS 

