/**
  * @file ComputeBasisFunctions_imp.hpp 
  *
  * @brief
  * implements functions to compute basis functions 
  *
  * @anchor _ComputeBasisFunctions_imp_hpp_  
  *
  * @author
  * Tanuharja R.A. @n 
  * contact: rezha.tanuharja@tum.de 
  */

#ifndef COMPUTE_BASIS_FUNCTIONS_IMPLEMENTATIONS 
#define COMPUTE_BASIS_FUNCTIONS_IMPLEMENTATIONS 

#ifndef BASIS_FUNCTIONS_DECLARATIONS 
  #include "BasisFunctions.hpp" 
#endif


namespace BasisFunctions {


template < typename Z > 
/**
  * @private 
  *
  * @brief  
  * Computes n! recursively. @n 
  * Used instead of boost implementation which return floating number. 
  *
  * @tparam Z non-negative integer e.g. size_t 
  */
Z Factorial ( const Z n ); 


template < typename Z, typename C >
/**
  * @private 
  *
  * @brief  
  * Overload multiplication operation between integer and complex number. @n
  *
  * @tparam Z non-negative integer e.g. size_t 
  * @tparam C complex float number e.g. std::complex<double>  
  */
C operator* ( const Z Integer, const C Complex );


template < typename Z, typename C > 
/**
  * @private 
  *
  * @brief  
  * Compute probabilist Hermite polynomial of given index and location.
  * The polynomial is NOT normalised. 
  *
  * @tparam Z non-negative integer e.g. size_t 
  * @tparam C complex float number e.g. std::complex<double>  
  *
  * @param Index indicates which polynomial to use 
  * @param X     argument to evaluate the polynomial 
  *
  * @return H_{Index} (X)
  */
C HermitePolynomial ( const Z Index, const C X );


} // namespace BasisFunctions 


namespace BasisFunctions {


template < typename Z, typename R, typename C > 
Vector<C> ComputeHermiteBasis (
  const Vector<Z>& Indices, const Vector<C>& Args, const Z SetSize 
) {


  #ifdef DEBUG_BASIS_FUNCTIONS 

  if ( SetSize < 1 || SetSize > 100 ) {
    throw std::runtime_error (
      "ComputeHermiteBasis: SetSize not in [1,100] range"
    );
  }

  #endif // DEBUG_BASIS_FUNCTIONS
  
  
  auto nProducts = Indices.size() / SetSize; 
  auto nSamples  = Args.size()    / SetSize; 


  #ifdef DEBUG_BASIS_FUNCTIONS 

  if ( SetSize * nProducts - Indices.size() != 0 ) {
    throw std::runtime_error (
      "ComputeHermiteBasis: Indices vector size not multiple of SetSize"
    );
  }

  if ( SetSize * nSamples - Args.size() != 0 ) {
    throw std::runtime_error (
      "ComputeHermiteBasis: Args vector size not multiple of SetSize"
    );
  }

  #endif // DEBUG_BASIS_FUNCTIONS 
  
  
  Vector<C> result;
  result.reserve ( nProducts * nSamples ); 

  for ( auto i = 0; i < nSamples; i++ ) {
    for ( auto j = 0; j < nProducts; j++ ) {

      C unity = 1.0; 

      result.push_back (
        std::transform_reduce (
          Indices.begin() + j * SetSize, 
          Indices.begin() + j * SetSize + SetSize, 
             Args.begin() + i * SetSize, 

          unity, 

          // Basis function is product of Hermite polynomials
          []( const auto N, const auto P ) { return N * P; }, 

          // P is normalised Hermite polynomial w/ given index and location
          []( const auto idx, const auto x ) {
            return HermitePolynomial<Z,C> ( idx, x ) /
              std::sqrt<R> ( Factorial<Z> ( idx )  );
          }

        )
      );

    }
  }

  return result; 

} // ComputeHermiteBasis 


} // namespace BasisFunctions 


namespace BasisFunctions {


template < typename Z > 
Z Factorial ( const Z n ) {


  #ifdef DEBUG_BASIS_FUNCTIONS 

  if ( n < 0 || n > 100 ) {
    throw std::runtime_error (
      "Factorial: n not in [0,100] range"
    );
  }

  #endif // DEBUG_BASIS_FUNCTIONS 
  
  
  if ( n == 0 ) { return 1; }

  return n * Factorial<Z> ( n - 1 ); 

} // Factorial 


template < typename Z, typename C >
C operator* ( const Z Integer, const C Complex ) {

  C result; 
  
  result.real( Integer * Complex.real() );
  result.imag( Integer * Complex.imag() );

  return result;

} // operator* 


template < typename Z, typename C > 
C HermitePolynomial ( const Z Index, const C X ) {


  #ifdef DEBUG_BASIS_FUNCTIONS 

  if ( Index < 0 || Index > 100 ) {
    throw std::runtime_error (
      "HermitePolynomial: Index not in [0,100] range"
    );
  }

  #endif // DEBUG_BASIS_FUNCTIONS 
   
  
  if ( Index == 0 ) { return 1.0; }
  if ( Index == 1 ) { return X;   }

  return (
    X             * HermitePolynomial<Z,C> ( Index - 1, X ) - 
    ( Index - 1 ) * HermitePolynomial<Z,C> ( Index - 2, X )
  );

} // HermitePolynomial 


} // namespace BasisFunctions 


#endif // COMPUTE_BASIS_FUNCTIONS_IMPLEMENTATIONS 

