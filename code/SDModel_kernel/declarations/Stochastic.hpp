/**
  * @file Stochastic.hpp 
  *
  * @brief
  * declares functions for stochastic model 
  *
  * @author
  * Tanuharja R.A. @n 
  * contact: rezha.tanuharja@tum.de 
  */

#ifndef STOCHASTIC_DECLARATIONS 
#define STOCHASTIC_DECLARATIONS 

#ifndef LIBRARIES_LOADER_RS 
  #include "LibrariesLoader_RS.hpp" 
#endif 

/**
  * @namespace Stochastic 
  * Contain functions for stochastic model 
  *
  * @anchor _Stochastic_ 
  */
namespace Stochastic {

typedef size_t          Z; 
typedef float           R; 
typedef std::complex<R> C; 

typedef std::vector<Z> VectorZ; 
typedef std::vector<R> VectorR; 
typedef std::vector<C> VectorC; 

typedef Eigen::Matrix<R, Eigen::Dynamic, Eigen::Dynamic> MatrixXR; 
typedef Eigen::Matrix<C, Eigen::Dynamic, Eigen::Dynamic> MatrixXC; 

typedef Eigen::Vector<R, Eigen::Dynamic> VectorXR; 
typedef Eigen::Vector<C, Eigen::Dynamic> VectorXC; 

/**
  * @namespace DirectMC  
  * Contain functions for Direct MC method 
  *
  * @anchor _DirectMC_ 
  */
namespace DirectMC {
} // namespace DirectMC 

} // namespace Stochastic 

#ifndef MASS_SPRING_DAMPER_DMC_DECLARATIONS 
  #include "MassSpringDamper_DMC.hpp" 
#endif 

#endif // STOCHASTIC_DECLARATIONS 

