/**
  * @file Normal_imp.hpp 
  *
  * @brief
  * implements template functions to generate RVs from normal distributions 
  *
  * @anchor _Normal_imp_hpp_ 
  *
  * @author
  * Tanuharja R.A. @n 
  * contact: rezha.tanuharja@tum.de 
  */

#ifndef NORMAL_IMPLEMENTATIONS 
#define NORMAL_IMPLEMENTATIONS 

#ifndef NORMAL_DECLARATIONS 
  #include "Normal.hpp" 
#endif 

namespace Statistics {

namespace Normal {


template < typename Z, typename R, typename C >
Vector<C> RandomSampling (
  const Z nPoints, const Z Dim, 
  const R Mean, const R StdDev
) {
  

  #ifdef DEBUG_STATISTICS 

  if ( nPoints < 1 || nPoints > 80000 ) {
    throw std::runtime_error (
      "Normal::RandomSampling: nPoints not in [1,80000] range"
    );
  }

  if ( Dim < 1 || Dim > 100 ) {
    throw std::runtime_error (
      "Normal::RandomSampling: Dim not in [1,100] range"
    );
  }

  if ( StdDev < 0 ) {
    throw std::runtime_error (
      "Normal::RandomSampling: negative StdDev"
    );
  }

  #endif // DEBUG_STATISTICS 
  

  std::random_device device; 
  std::mt19937 generator ( device () );

  std::normal_distribution<R> RandomVariable ( Mean, StdDev );

  Vector<C> result ( nPoints * Dim ); 

  for ( auto i = 0; i < result.size(); i++ ) {
    result[i] = C ( RandomVariable ( generator ) );
  }

  return result;

} // RandomSampling 


} // namespace Normal 

} // namespace Statistics 

#endif //NORMAL_IMPLEMENTATIONS 

