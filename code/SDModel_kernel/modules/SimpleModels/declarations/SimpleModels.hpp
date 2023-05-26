/**
  * @file SimpleModels.hpp 
  *
  * @brief
  * declares template functions to solve simple models 
  *
  * @author
  * Tanuharja R.A. @n 
  * contact: rezha.tanuharja@tum.de 
  */

#ifndef SIMPLE_MODELS_DECLARATIONS 
#define SIMPLE_MODELS_DECLARATIONS 

#ifndef LIBRARIES_LOADER_SM 
  #include "LibrariesLoader_SM.hpp" 
#endif 

/**
  * @namespace SimpleModels 
  * Contain functions to solve simple models 
  *
  * @anchor _SimpleModels_ 
  */
namespace SimpleModels {

template < typename T > 
using Vector = std::vector<T>; 

template < typename T > 
using MatrixXT = Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic >; 

template < typename T > 
using VectorXT = Eigen::Vector< T, Eigen::Dynamic >; 

} // namespace SimpleModels 

#ifndef MASS_SPRING_DAMPER_DECLARATIONS 
  #include "MassSpringDamper.hpp" 
#endif 

#endif // SIMPLE_MODELS_DECLARATIONS 

