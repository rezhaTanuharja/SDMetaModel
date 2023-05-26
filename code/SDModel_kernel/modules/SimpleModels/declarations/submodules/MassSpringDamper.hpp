/**
  * @file MassSpringDamper.hpp 
  *
  * @brief
  * declares template functions to solve mass-spring-damper models 
  *
  * @author
  * Tanuharja R.A. @n 
  * contact: rezha.tanuharja@tum.de 
  */

#ifndef MASS_SPRING_DAMPER_DECLARATIONS 
#define MASS_SPRING_DAMPER_DECLARATIONS 

#ifndef SIMPLE_MODELS_DECLARATIONS 
  #include "SimpleModels.hpp" 
#endif 

namespace SimpleModels {

/**
  * @namespace MassSpringDamper 
  * Contain functions to solve mass-spring-damper models 
  *
  * @anchor _MassSpringDamper_ 
  */
namespace MassSpringDamper {


template < typename R > 
/**
  * @brief 
  * Compute mass matrix for given point masses. @n 
  * Implemented in @ref _MassSpringDamper_imp_hpp_ 
  *
  * @tparam R float number e.g. double 
  *
  * @param Masses vector of point masses (size multiples of Dim)
  *
  * @return mass matrix in column major storage order 
  */
Vector<R> MassMatrix ( const Vector<R>& Masses );


template < typename R > 
/**
  * @brief 
  * Compute stiffness matrix for given spring constants. @n 
  * Implemented in @ref _MassSpringDamper_imp_hpp_ 
  *
  * @tparam R float number e.g. double 
  *
  * @param Springs vector of spring constants (size multiples of Dim)
  *
  * @return stiffness matrix in column major storage order 
  */
Vector<R> StiffnessMatrix ( const Vector<R>& Springs );


template < typename R, typename C > 
/**
  * @brief 
  * Compute dynamic stiffness matrix of mass-spring-damper model. @n
  * Implemented in @ref _MassSpringDamper_imp_hpp_ 
  *
  * @tparam R float number e.g. double 
  * @tparam C complex float number e.g. std::complex<double>  
  *
  * @param MassMatrix      mass      matrix in row/column major storage order 
  * @param StiffnessMatrix stiffness matrix in row/column major storage order 
  * @param DampingRatio    structural parameter 
  * @param Omega           angular speed 
  *
  * @return K - Omega** M + jOmega C 
  */
Vector<C> DynamicStiffnessMatrix (
  const Vector<R>& MassMatrix, 
  const Vector<R>& StiffnessMatrix, 
  const R DampingRatio, 
  const R Omega 
);


} // namespace MassSpringDamper 

} // namespace SimpleModels 


#ifndef MASS_SPRING_DAMPER_IMPLEMENTATIONS 
  #include "MassSpringDamper_imp.hpp" 
#endif 

#endif // MASS_SPRING_DAMPER_DECLARATIONS 

