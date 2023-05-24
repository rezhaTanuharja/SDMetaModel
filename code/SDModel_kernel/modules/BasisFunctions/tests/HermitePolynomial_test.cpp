/**
  * @file HermitePolynomial_test.cpp 
  *
  * @brief
  * test functions to compute Hermite polynomials, 
  * including errors thrown in DEBUG mode
  *
  * @author
  * Tanuharja R.A. @n 
  * contact: rezha.tanuharja@tum.de 
  */

#include "BasisFunctions.hpp" 
#include <gtest/gtest.h> 


TEST ( HermitePolynomial, H_0_ComplexInput ) {

  typedef std::complex<double> Complex; 

  size_t index = 0; 

  Complex X        ( 3.0, 2.0 ); 
  Complex expected ( 1.0, 0.0 ); 

  auto result = BasisFunctions::HermitePolynomial ( index, X );

  EXPECT_FLOAT_EQ ( result.real(), expected.real() );
  EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

} // HermiteBasis, H_0_ComplexInput 


TEST ( HermitePolynomial, H_0_RealInput ) {

  typedef std::complex<double> Complex; 

  size_t index = 0; 

  Complex X        ( 1.0, 0.0 ); 
  Complex expected ( 1.0, 0.0 ); 

  auto result = BasisFunctions::HermitePolynomial ( index, X );

  EXPECT_FLOAT_EQ ( result.real(), expected.real() );
  EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

} // HermiteBasis, H_0_RealInput  


TEST ( HermitePolynomial, H_0_ImagInput ) {

  typedef std::complex<double> Complex; 

  size_t index = 0; 

  Complex X        ( 0.0, 5.0 ); 
  Complex expected ( 1.0, 0.0 ); 

  auto result = BasisFunctions::HermitePolynomial ( index, X );

  EXPECT_FLOAT_EQ ( result.real(), expected.real() );
  EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

} // HermiteBasis, H_0_ImagInput  


TEST ( HermitePolynomial, H_0_NegInput ) {

  typedef std::complex<double> Complex; 

  size_t index = 0; 

  Complex X        (-7.0,10.0 ); 
  Complex expected ( 1.0, 0.0 ); 

  auto result = BasisFunctions::HermitePolynomial ( index, X );

  EXPECT_FLOAT_EQ ( result.real(), expected.real() );
  EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

} // HermiteBasis, H_0_NegInput 


TEST ( HermitePolynomial, H_1_ComplexInput ) {

  typedef std::complex<double> Complex; 

  size_t index = 1; 

  Complex X        ( 3.0, 2.0 ); 
  Complex expected = X; 

  auto result = BasisFunctions::HermitePolynomial ( index, X );

  EXPECT_FLOAT_EQ ( result.real(), expected.real() );
  EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

} // HermiteBasis, H_1_ComplexInput 


TEST ( HermitePolynomial, H_1_RealInput ) {

  typedef std::complex<double> Complex; 

  size_t index = 1; 

  Complex X        ( 1.0, 0.0 ); 
  Complex expected = X; 

  auto result = BasisFunctions::HermitePolynomial ( index, X );

  EXPECT_FLOAT_EQ ( result.real(), expected.real() );
  EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

} // HermiteBasis, H_1_RealInput  


TEST ( HermitePolynomial, H_1_ImagInput ) {

  typedef std::complex<double> Complex; 

  size_t index = 1; 

  Complex X        ( 0.0, 5.0 ); 
  Complex expected = X; 

  auto result = BasisFunctions::HermitePolynomial ( index, X );

  EXPECT_FLOAT_EQ ( result.real(), expected.real() );
  EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

} // HermiteBasis, H_1_ImagInput  


TEST ( HermitePolynomial, H_1_NegInput ) {

  typedef std::complex<double> Complex; 

  size_t index = 1; 

  Complex X        (-7.0,10.0 ); 
  Complex expected = X; 

  auto result = BasisFunctions::HermitePolynomial ( index, X );

  EXPECT_FLOAT_EQ ( result.real(), expected.real() );
  EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

} // HermiteBasis, H_1_NegInput 


TEST ( HermitePolynomial, H_2_ComplexInput ) {

  typedef std::complex<double> Complex; 

  size_t index = 2; 

  Complex X        ( 3.0, 2.0 ); 
  Complex expected = X * X - 1.0; 

  auto result = BasisFunctions::HermitePolynomial ( index, X );

  EXPECT_FLOAT_EQ ( result.real(), expected.real() );
  EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

} // HermiteBasis, H_2_ComplexInput 


TEST ( HermitePolynomial, H_2_RealInput ) {

  typedef std::complex<double> Complex; 

  size_t index = 2; 

  Complex X        ( 1.0, 0.0 ); 
  Complex expected = X * X - 1.0; 

  auto result = BasisFunctions::HermitePolynomial ( index, X );

  EXPECT_FLOAT_EQ ( result.real(), expected.real() );
  EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

} // HermiteBasis, H_2_RealInput  


TEST ( HermitePolynomial, H_2_ImagInput ) {

  typedef std::complex<double> Complex; 

  size_t index = 2; 

  Complex X        ( 0.0, 5.0 ); 
  Complex expected = X * X - 1.0; 

  auto result = BasisFunctions::HermitePolynomial ( index, X );

  EXPECT_FLOAT_EQ ( result.real(), expected.real() );
  EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

} // HermiteBasis, H_2_ImagInput  


TEST ( HermitePolynomial, H_2_NegInput ) {

  typedef std::complex<double> Complex; 

  size_t index = 2; 

  Complex X        (-7.0,10.0 ); 
  Complex expected = X * X - 1.0; 

  auto result = BasisFunctions::HermitePolynomial ( index, X );

  EXPECT_FLOAT_EQ ( result.real(), expected.real() );
  EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

} // HermiteBasis, H_2_NegInput 


TEST ( HermitePolynomial, H_3_ComplexInput ) {

  typedef std::complex<double> Complex; 

  size_t index = 3; 

  Complex X        ( 3.0, 2.0 ); 
  Complex expected = X * X * X - 3.0 * X; 

  auto result = BasisFunctions::HermitePolynomial ( index, X );

  EXPECT_FLOAT_EQ ( result.real(), expected.real() );
  EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

} // HermiteBasis, H_3_ComplexInput 


TEST ( HermitePolynomial, H_3_RealInput ) {

  typedef std::complex<double> Complex; 

  size_t index = 3; 

  Complex X        ( 1.0, 0.0 ); 
  Complex expected = X * X * X - 3.0 * X; 

  auto result = BasisFunctions::HermitePolynomial ( index, X );

  EXPECT_FLOAT_EQ ( result.real(), expected.real() );
  EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

} // HermiteBasis, H_3_RealInput  


TEST ( HermitePolynomial, H_3_ImagInput ) {

  typedef std::complex<double> Complex; 

  size_t index = 3; 

  Complex X        ( 0.0, 5.0 ); 
  Complex expected = X * X * X - 3.0 * X; 

  auto result = BasisFunctions::HermitePolynomial ( index, X );

  EXPECT_FLOAT_EQ ( result.real(), expected.real() );
  EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

} // HermiteBasis, H_3_ImagInput  


TEST ( HermitePolynomial, H_3_NegInput ) {

  typedef std::complex<double> Complex; 

  size_t index = 3; 

  Complex X        (-7.0,10.0 ); 
  Complex expected = X * X * X - 3.0 * X; 

  auto result = BasisFunctions::HermitePolynomial ( index, X );

  EXPECT_FLOAT_EQ ( result.real(), expected.real() );
  EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

} // HermiteBasis, H_3_NegInput 


TEST ( HermitePolynomial, H_4_ComplexInput ) {

  typedef std::complex<double> Complex; 

  size_t index = 4; 

  Complex X        ( 3.0, 2.0 ); 
  Complex expected = X * X * X * X - 6.0 * X * X + 3.0; 

  auto result = BasisFunctions::HermitePolynomial ( index, X );

  EXPECT_FLOAT_EQ ( result.real(), expected.real() );
  EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

} // HermiteBasis, H_4_ComplexInput 


TEST ( HermitePolynomial, H_4_RealInput ) {

  typedef std::complex<double> Complex; 

  size_t index = 4; 

  Complex X        ( 1.0, 0.0 ); 
  Complex expected = X * X * X * X - 6.0 * X * X + 3.0; 

  auto result = BasisFunctions::HermitePolynomial ( index, X );

  EXPECT_FLOAT_EQ ( result.real(), expected.real() );
  EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

} // HermiteBasis, H_4_RealInput  


TEST ( HermitePolynomial, H_4_ImagInput ) {

  typedef std::complex<double> Complex; 

  size_t index = 4; 

  Complex X        ( 0.0, 5.0 ); 
  Complex expected = X * X * X * X - 6.0 * X * X + 3.0; 

  auto result = BasisFunctions::HermitePolynomial ( index, X );

  EXPECT_FLOAT_EQ ( result.real(), expected.real() );
  EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

} // HermiteBasis, H_4_ImagInput  


TEST ( HermitePolynomial, H_4_NegInput ) {

  typedef std::complex<double> Complex; 

  size_t index = 4; 

  Complex X        (-7.0,10.0 ); 
  Complex expected = X * X * X * X - 6.0 * X * X + 3.0; 

  auto result = BasisFunctions::HermitePolynomial ( index, X );

  EXPECT_FLOAT_EQ ( result.real(), expected.real() );
  EXPECT_FLOAT_EQ ( result.imag(), expected.imag() );

} // HermiteBasis, H_4_NegInput 


#ifdef DEBUG_BASIS_FUNCTIONS 


TEST ( HermitePolynomial, H_10_NegIndex ) {

  typedef std::complex<double> Complex; 

  size_t index = -1; 

  Complex X        (-7.0,10.0 ); 

  bool exception_thrown = false; 

  try {

    [[maybe_unused]] auto result = 
      BasisFunctions::HermitePolynomial ( index, X );

  } catch ( const std::exception& e ) {

    EXPECT_STREQ (
      
      "HermitePolynomial: Index not in [0,100] range", e.what()
      
    );

    exception_thrown = true;

  }

  EXPECT_TRUE ( exception_thrown );

} // HermiteBasis, H_10_NegIndex 


TEST ( HermitePolynomial, H_10_HugeIndex ) {

  typedef std::complex<double> Complex; 

  size_t index = 1000; 

  Complex X        (-7.0,10.0 ); 

  bool exception_thrown = false; 

  try {

    [[maybe_unused]] auto result = 
      BasisFunctions::HermitePolynomial ( index, X );

  } catch ( const std::exception& e ) {

    EXPECT_STREQ (
      
      "HermitePolynomial: Index not in [0,100] range", e.what()
      
    );

    exception_thrown = true;

  }

  EXPECT_TRUE ( exception_thrown );

} // HermiteBasis, H_10_HugeIndex 


#endif // DEBUG_BASIS_FUNCTIONS 

