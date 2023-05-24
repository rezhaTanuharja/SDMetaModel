/**
  * @file HermiteBasis_test.cpp 
  *
  * @brief
  * test functions to compute basis from products of Hermite polynomials, 
  * including errors thrown in DEBUG mode
  *
  * @author
  * Tanuharja R.A. @n 
  * contact: rezha.tanuharja@tum.de 
  */

#include "BasisFunctions.hpp" 
#include <gtest/gtest.h> 


TEST ( HermiteBasis, H_0000_4_ComplexInput ) {

  typedef std::complex<double> Complex;
  typedef std::vector<Complex> Vector;

  std::vector<size_t> indices { 0, 0, 0, 0 };

  Vector X {

    Complex (  2.0,  3.0 ), 
    Complex ( -5.0,  2.0 ), 
    Complex (  0.0,  2.0 ), 
    Complex ( -5.0,  0.0 ) 

  };

  size_t dim = 4;

  auto result =
    BasisFunctions::ComputeHermiteBasis<
      
      size_t, double, std::complex<double>
      
    > ( indices, X, dim );

  Complex expected;
  expected.real(  1.0 );
  expected.imag(  0.0 );

  ASSERT_EQ ( result.size(), 1 );
  EXPECT_FLOAT_EQ ( result[0].real(), expected.real() );
  EXPECT_FLOAT_EQ ( result[0].imag(), expected.imag() );

} // HermiteBasis, H_0000_4_ComplexInput 


TEST ( HermiteBasis, H_2000_4_ComplexInput ) {

  typedef std::complex<double> Complex;
  typedef std::vector<Complex> Vector;

  std::vector<size_t> indices { 2, 0, 0, 0 };

  Vector X {

    Complex (  2.0,  3.0 ), 
    Complex ( -5.0,  2.0 ), 
    Complex (  0.0,  2.0 ), 
    Complex ( -5.0,  0.0 ) 

  };

  size_t dim = 4;

  auto result =
    BasisFunctions::ComputeHermiteBasis<
      
      size_t, double, std::complex<double>
      
    > ( indices, X, dim );

  Complex expected = ( X[0] * X[0] - 1.0 ) / std::sqrt(2);

  ASSERT_EQ ( result.size(), 1 );
  EXPECT_FLOAT_EQ ( result[0].real(), expected.real() );
  EXPECT_FLOAT_EQ ( result[0].imag(), expected.imag() );

} // HermiteBasis, H_2000_4_ComplexInput 


TEST ( HermiteBasis, H_1010_4_ComplexInput ) {

  typedef std::complex<double> Complex;
  typedef std::vector<Complex> Vector;

  std::vector<size_t> indices { 1, 0, 1, 0 };

  Vector X {

    Complex (  2.0,  3.0 ), 
    Complex ( -5.0,  2.0 ), 
    Complex (  0.0,  2.0 ), 
    Complex ( -5.0,  0.0 ) 

  };

  size_t dim = 4;

  auto result =
    BasisFunctions::ComputeHermiteBasis<
      
      size_t, double, std::complex<double>
      
    > ( indices, X, dim );

  Complex expected = X[0] * X[2];

  ASSERT_EQ ( result.size(), 1 );
  EXPECT_FLOAT_EQ ( result[0].real(), expected.real() );
  EXPECT_FLOAT_EQ ( result[0].imag(), expected.imag() );

} // HermiteBasis, H_1010_4_ComplexInput 


TEST ( HermiteBasis, H_0003_4_ComplexInput ) {

  typedef std::complex<double> Complex;
  typedef std::vector<Complex> Vector;

  std::vector<size_t> indices { 0, 0, 0, 3 };

  Vector X {

    Complex (  2.0,  3.0 ), 
    Complex ( -5.0,  2.0 ), 
    Complex (  0.0,  2.0 ), 
    Complex ( -5.0,  0.0 ) 

  };

  size_t dim = 4;

  auto result =
    BasisFunctions::ComputeHermiteBasis<
      
      size_t, double, std::complex<double>
      
    > ( indices, X, dim );

  Complex expected = ( X[3] * X[3] * X[3] - X[3] * 3.0 ) / std::sqrt(6);

  ASSERT_EQ ( result.size(), 1 );
  EXPECT_FLOAT_EQ ( result[0].real(), expected.real() );
  EXPECT_FLOAT_EQ ( result[0].imag(), expected.imag() );

} // HermiteBasis, H_0003_4_ComplexInput 


TEST ( HermiteBasis, H_211203_2_ComplexInput ) {

  typedef std::complex<double> Complex;
  typedef std::vector<Complex> Vector;

  std::vector<size_t> indices { 2, 1, 1, 2, 0, 3 };

  Vector X {

    Complex (  2.0,  3.0 ), 
    Complex ( -5.0,  2.0 ) 

  };

  size_t dim = 2;

  auto result =
    BasisFunctions::ComputeHermiteBasis<
      
      size_t, double, std::complex<double>
      
    > ( indices, X, dim );

  std::vector<std::complex<double>> expected;

  expected.push_back (
    X[1] * ( X[0] * X[0] - 1.0 ) / std::sqrt(2)
  );

  expected.push_back (
    X[0] * ( X[1] * X[1] - 1.0 ) / std::sqrt(2)
  );

  expected.push_back (
    ( X[1] * X[1] * X[1] - X[1] * 3.0 ) / std::sqrt(6)
  );

  ASSERT_EQ ( result.size(), expected.size() );

  for ( auto i = 0; i < expected.size(); i++ ) {
    
    EXPECT_FLOAT_EQ ( result[i].real(), expected[i].real() );
    EXPECT_FLOAT_EQ ( result[i].imag(), expected[i].imag() );

  }

} // HermiteBasis, H_211203_2_ComplexInput 


TEST ( HermiteBasis, H_221213_2_ComplexInput ) {

  typedef std::complex<double> Complex;
  typedef std::vector<Complex> Vector;

  std::vector<size_t> indices { 2, 2, 1, 2, 1, 3 };

  Vector X {

    Complex (  2.0,  3.0 ), 
    Complex ( -5.0,  2.0 ), 
    Complex (  7.0,  7.0 ), 
    Complex (  5.0, -4.0 ) 

  };

  size_t dim = 2;

  auto result =
    BasisFunctions::ComputeHermiteBasis<
      
      size_t, double, std::complex<double>
      
    > ( indices, X, dim );

  std::vector<std::complex<double>> expected;

  expected.push_back (
    ( X[0] * X[0] - 1.0 ) * ( X[1] * X[1] - 1.0 ) / 2.0
  );

  expected.push_back (
    X[0] * ( X[1] * X[1] - 1.0 ) / std::sqrt(2)
  );

  expected.push_back (
    X[0] * ( X[1] * X[1] * X[1] - X[1] * 3.0 ) / std::sqrt(6)
  );

  expected.push_back (
    ( X[2] * X[2] - 1.0 ) * ( X[3] * X[3] - 1.0 ) / 2.0
  );

  expected.push_back (
    X[2] * ( X[3] * X[3] - 1.0 ) / std::sqrt(2)
  );

  expected.push_back (
    X[2] * ( X[3] * X[3] * X[3] - X[3] * 3.0 ) / std::sqrt(6)
  );

  ASSERT_EQ ( result.size(), expected.size() );

  for ( auto i = 0; i < expected.size(); i++ ) {
    
    EXPECT_FLOAT_EQ ( result[i].real(), expected[i].real() );
    EXPECT_FLOAT_EQ ( result[i].imag(), expected[i].imag() );

  }

} // HermiteBasis, H_221213_2_ComplexInput 


#ifdef DEBUG_BASIS_FUNCTIONS 


TEST ( HermiteBasis, H_5312_4_ZeroSetSize ) {

  typedef std::complex<double> Complex;
  typedef std::vector<Complex> Vector;

  std::vector<size_t> indices { 5, 3, 1, 2 };

  Vector X {

    Complex (  2.0,  3.0 ), 
    Complex ( -5.0,  2.0 ), 
    Complex (  0.0,  2.0 ), 
    Complex ( -5.0,  0.0 ) 

  };

  size_t dim = 0;

  bool exception_thrown = false; 

  try {
    
    auto result =
      BasisFunctions::ComputeHermiteBasis<
        
        size_t, double, std::complex<double>
        
      > ( indices, X, dim );
    
  } catch ( const std::exception& e ) {
    
    EXPECT_STREQ (
      
      "ComputeHermiteBasis: SetSize not in [1,100] range", e.what()
      
    );

    exception_thrown = true;
    
  }

  EXPECT_TRUE ( exception_thrown );

} // HermiteBasis, H_5312_4_ZeroSetSize


TEST ( HermiteBasis, H_5312_4_NegSetSize ) {

  typedef std::complex<double> Complex;
  typedef std::vector<Complex> Vector;

  std::vector<size_t> indices { 5, 3, 1, 2 };

  Vector X {

    Complex (  2.0,  3.0 ), 
    Complex ( -5.0,  2.0 ), 
    Complex (  0.0,  2.0 ), 
    Complex ( -5.0,  0.0 ) 

  };

  size_t dim = -1;

  bool exception_thrown = false; 

  try {
    
    auto result =
      BasisFunctions::ComputeHermiteBasis<
        
        size_t, double, std::complex<double>
        
      > ( indices, X, dim );
    
  } catch ( const std::exception& e ) {
    
    EXPECT_STREQ (
      
      "ComputeHermiteBasis: SetSize not in [1,100] range", e.what()
      
    );

    exception_thrown = true;
    
  }

  EXPECT_TRUE ( exception_thrown );

} // HermiteBasis, H_5312_4_NegSetSize


TEST ( HermiteBasis, H_5312_4_HugeSetSize ) {

  typedef std::complex<double> Complex;
  typedef std::vector<Complex> Vector;

  std::vector<size_t> indices { 5, 3, 1, 2 };

  Vector X {

    Complex (  2.0,  3.0 ), 
    Complex ( -5.0,  2.0 ), 
    Complex (  0.0,  2.0 ), 
    Complex ( -5.0,  0.0 ) 

  };

  size_t dim = 1000;

  bool exception_thrown = false; 

  try {
    
    auto result =
      BasisFunctions::ComputeHermiteBasis<
        
        size_t, double, std::complex<double>
        
      > ( indices, X, dim );
    
  } catch ( const std::exception& e ) {
    
    EXPECT_STREQ (
      
      "ComputeHermiteBasis: SetSize not in [1,100] range", e.what()
      
    );

    exception_thrown = true;
    
  }

  EXPECT_TRUE ( exception_thrown );

} // HermiteBasis, H_5312_4_HugeSetSize


TEST ( HermiteBasis, H_33123_2_WrongIndicesSize ) {

  typedef std::complex<double> Complex;
  typedef std::vector<Complex> Vector;

  std::vector<size_t> indices { 3, 3, 1, 2, 3 };

  Vector X {

    Complex (  2.0,  3.0 ), 
    Complex ( -5.0,  2.0 ), 
    Complex (  0.0,  2.0 ), 
    Complex ( -5.0,  0.0 ) 

  };

  size_t dim = 2;

  bool exception_thrown = false; 

  try {
    
    auto result =
      BasisFunctions::ComputeHermiteBasis<
        
        size_t, double, std::complex<double>
        
      > ( indices, X, dim );
    
  } catch ( const std::exception& e ) {
    
    EXPECT_STREQ (
      
      "ComputeHermiteBasis: Indices vector size not multiple of SetSize", 
      e.what()
      
    );

    exception_thrown = true;
    
  }

  EXPECT_TRUE ( exception_thrown );

} // HermiteBasis, H_33123_2_HugeSetSize


TEST ( HermiteBasis, H_3312_2_WrongArgsSize ) {

  typedef std::complex<double> Complex;
  typedef std::vector<Complex> Vector;

  std::vector<size_t> indices { 3, 3, 1, 2 };

  Vector X {

    Complex (  2.0,  3.0 ), 
    Complex ( -5.0,  2.0 ), 
    Complex (  0.0,  2.0 ), 
    Complex ( -5.0,  0.0 ), 
    Complex ( -5.0,  0.0 ) 

  };

  size_t dim = 2;

  bool exception_thrown = false; 

  try {
    
    auto result =
      BasisFunctions::ComputeHermiteBasis<
        
        size_t, double, std::complex<double>
        
      > ( indices, X, dim );
    
  } catch ( const std::exception& e ) {
    
    EXPECT_STREQ (
      
      "ComputeHermiteBasis: Args vector size not multiple of SetSize", 
      e.what()
      
    );

    exception_thrown = true;
    
  }

  EXPECT_TRUE ( exception_thrown );

} // HermiteBasis, H_33123_2_HugeSetSize


#endif // DEBUG_BASIS_FUNCTIONS 

