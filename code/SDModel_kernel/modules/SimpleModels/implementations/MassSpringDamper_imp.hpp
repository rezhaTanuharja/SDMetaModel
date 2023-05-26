/**
  * @file MassSpringDamper_imp.hpp 
  *
  * @brief
  * implements template functions to solve mass-spring-damper models 
  *
  * @anchor _MassSpringDamper_imp_hpp_ 
  *
  * @author
  * Tanuharja R.A. @n 
  * contact: rezha.tanuharja@tum.de 
  */

#ifndef MASS_SPRING_DAMPER_IMPLEMENTATIONS 
#define MASS_SPRING_DAMPER_IMPLEMENTATIONS 

#ifndef MASS_SPRING_DAMPER_DECLARATIONS 
  #include "MassSpringDamper.hpp" 
#endif 

namespace SimpleModels::MassSpringDamper {


template < typename R > 
Vector<R> MassMatrix ( const Vector<R>& Masses ) {

  auto Dim = Masses.size(); 

  Vector<R> result ( Dim * Dim, 0.0 ); 

  for ( auto i = 0; i < Dim; i++ ) {
    result[i+i*Dim] = Masses[i];
  }

  return result; 

} // MassMatrix 


template < typename R > 
Vector<R> StiffnessMatrix ( const Vector<R>& Springs ) {

  auto Dim = Springs.size(); 

  if ( Dim == 0 ) {
    return Vector<R> (Springs); 
  }

  Vector<R> result ( Dim * Dim, 0.0 ); 

  result[0] = Springs[0] + Springs[1]; 
  result[1] =            - Springs[1]; 

  for ( auto i = 1; i < Dim - 1; i++ ) {

    result[i+i*Dim-1] = -Springs[i]; 
    result[i+i*Dim  ] =  Springs[i] + Springs[i+1]; 
    result[i+i*Dim+1] =             - Springs[i+1]; 

  }

  result[Dim*Dim-2] = -Springs[Dim-1]; 
  result[Dim*Dim-1] =  Springs[Dim-1]; 

  return result; 

} // StiffnessMatrix 


template < typename R, typename C > 
Vector<C> DynamicStiffnessMatrix (
  const Vector<R>& MassMatrix, 
  const Vector<R>& StiffnessMatrix, 
  const R DampingRatio, 
  const R Omega 
) {


  #ifdef DEBUG_SIMPLE_MODELS 

  if ( DampingRatio < 0 ) {
    throw std::runtime_error (
      "MassSpringDamper::DynamicStiffnessMatrix: negative DampingRatio"
    );
  }

  if ( MassMatrix.size() != StiffnessMatrix.size() ) {
    throw std::runtime_error (
      "MassSpringDamper::DynamicStiffnessMatrix: matrices size doesn't match"
    );
  }

  #endif // DEBUG_SIMPLE_MODELS 
  
  
  Vector<C> result ( MassMatrix.size() ); 

  std::transform (

    MassMatrix.begin(), 
    MassMatrix.end(), 

    StiffnessMatrix.begin(), 

    result.begin(), 

    [DampingRatio,Omega]( const auto m, const auto k ) {
      return C ( k - Omega*Omega*m, 2.0 * Omega * DampingRatio * k ); 
    }

  );
  
  return result; 

} // DynamicStiffnessMatrix 


} // namespace SimpleModels::MassSpringDamper 

#endif // MASS_SPRING_DAMPER_IMPLEMENTATIONS 

