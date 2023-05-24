/**
  * @file BasisFunctions.hpp 
  *
  * @brief
  * declares template functions to generate Hermite basis functions 
  *
  * @author
  * Tanuharja R.A. @n 
  * contact: rezha.tanuharja@tum.de 
  */

#ifndef BASIS_FUNCTIONS_DECLARATIONS 
#define BASIS_FUNCTIONS_DECLARATIONS 

#ifndef LIBRARIES_LOADER_BF 
  #include "LibrariesLoader_BF.hpp" 
#endif 

/** 
  * @namespace BasisFunctions 
  * 
  * @brief 
  * Contains functions to generate Hermite basis functions 
  * 
  * @anchor _BasisFunctions_ 
  */
namespace BasisFunctions {

template < typename T > 
using Vector = std::vector<T>; 


template < typename Z > 
/**
  * @brief 
  * Generate sets of indices of Hermite polynomials. 
  * Sorted with increasing sum of indices in each set. @n 
  * Implemented in @ref _GenerateIndices_imp_hpp_ 
  *
  * @tparam Z non-negative integer e.g. size_t 
  *
  * @param SetSize number of indices in each set 
  * @param MaxSum  largest allowable sum of indices in each set 
  *
  * @return vector of indices, size of which is a multiple of SetSize 
  */
Vector<Z> GenerateIndices ( const Z SetSize, const Z MaxSum ); 


template < typename Z > 
/**
  * @brief 
  * Remove sets of indices that contain large index. 
  * Operations are performed in-place, no new vector is generated. @n
  * Implemented in @ref _GenerateIndices_imp_hpp_ 
  *
  * @tparam Z non-negative integer e.g. size_t 
  *
  * @param Indices vector of indices e.g. from function GenerateIndices 
  * @param SetSize number of indices in each set 
  * @param iMax    largest allowable index 
  */
void RemoveLargeIndices ( Vector<Z>& Indices, const Z SetSize, const Z iMax );


template < typename Z, typename R, typename C > 
/**
  * @brief 
  * Compute basis functions from products of Hermite polynomials. @n 
  * Implemented in @ref _ComputeBasisFunctions_imp_hpp_ 
  *
  * @tparam Z non-negative integer e.g. size_t 
  * @tparam R float number e.g. double 
  * @tparam C complex float number e.g. std::complex<double>  
  *
  * @param Indices vector of indices e.g. from function GenerateIndices 
  * @param Args    vector of arguments for Hermite polynomials 
  * @param SetSize number of indices in each set 
  *
  * @return vector of products of Hermite polynomials 
  */
Vector<C> ComputeHermiteBasis (
  const Vector<Z>& Indices, const Vector<C>& Args, const Z SetSize 
);


} // namespace BasisFunctions 

#ifndef GENERATE_INDICES_IMPLEMENTATIONS 
  #include "GenerateIndices_imp.hpp"
#endif 

#ifndef COMPUTE_BASIS_FUNCTIONS_IMPLEMENTATIONS 
  #include "ComputeBasisFunctions_imp.hpp"
#endif 

#endif // BASIS_FUNCTIONS_DECLARATIONS 

