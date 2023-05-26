/**
  * @file DynamicStiffnessMatrix_test.cpp 
  *
  * @brief
  * test functions to compute dyn stiffness for mass-spring-damper models 
  * including errors thrown in DEBUG mode
  *
  * @author
  * Tanuharja R.A. @n 
  * contact: rezha.tanuharja@tum.de 
  */

#include "SimpleModels.hpp" 
#include <gtest/gtest.h> 


TEST ( DynamicStiffnessMatrix, NDOF_1_NoDamping ) {
  
  typedef double               R; 
  typedef std::complex<double> C; 
  
  std::vector<R> MassMatrix {
     15.0
  };
  
  std::vector<R> StiffnessMatrix {
    1000.0
  };
  
  R DampingRatio = 0.0; 
  R Omega        = 5.0; 
  
  std::vector<C> expected {
    C ( 625.0, 0.0 )
  };
  
  auto result = 
    SimpleModels::MassSpringDamper::DynamicStiffnessMatrix<R,C> (
      MassMatrix, StiffnessMatrix, 
      DampingRatio, Omega
    );
  
  ASSERT_EQ ( expected.size(), result.size() );

  for ( auto i = 0; i < result.size(); i++ ) {
    EXPECT_EQ ( expected[i].real(), result[i].real() );
    EXPECT_EQ ( expected[i].imag(), result[i].imag() );
  }

} // DynamicStiffnessMatrix, NDOF_1_NoDamping 


TEST ( DynamicStiffnessMatrix, NDOF_1_WithDamping ) {
  
  typedef double               R; 
  typedef std::complex<double> C; 
  
  std::vector<R> MassMatrix {
     15.0
  };
  
  std::vector<R> StiffnessMatrix {
    1000.0
  };
  
  R DampingRatio = 0.1; 
  R Omega        = 5.0; 
  
  std::vector<C> expected {
    C ( 625.0, 1000.0 )
  };
  
  auto result = 
    SimpleModels::MassSpringDamper::DynamicStiffnessMatrix<R,C> (
      MassMatrix, StiffnessMatrix, 
      DampingRatio, Omega
    );
  
  ASSERT_EQ ( expected.size(), result.size() );

  for ( auto i = 0; i < result.size(); i++ ) {
    EXPECT_EQ ( expected[i].real(), result[i].real() );
    EXPECT_EQ ( expected[i].imag(), result[i].imag() );
  }

} // DynamicStiffnessMatrix, NDOF_1_WithDamping 


TEST ( DynamicStiffnessMatrix, NDOF_2_NoDamping ) {
  
  typedef double               R; 
  typedef std::complex<double> C; 
  
  std::vector<R> MassMatrix {
    10.0,  0.0,
     0.0, 15.0
  };
  
  std::vector<R> StiffnessMatrix {
    1500.0, -500.0, 
    -500.0,  500.0 
  };
  
  R DampingRatio = 0.0; 
  R Omega        = 5.0; 
  
  std::vector<C> expected {
    C ( 1250.0, 0.0 ), C ( -500.0, 0.0 ), 
    C ( -500.0, 0.0 ), C (  125.0, 0.0 ) 
  };
  
  auto result = 
    SimpleModels::MassSpringDamper::DynamicStiffnessMatrix<R,C> (
      MassMatrix, StiffnessMatrix, 
      DampingRatio, Omega
    );
  
  ASSERT_EQ ( expected.size(), result.size() );

  for ( auto i = 0; i < result.size(); i++ ) {
    EXPECT_EQ ( expected[i].real(), result[i].real() );
    EXPECT_EQ ( expected[i].imag(), result[i].imag() );
  }

} // DynamicStiffnessMatrix, NDOF_2_NoDamping 


TEST ( DynamicStiffnessMatrix, NDOF_2_WithDamping ) {
  
  typedef double               R; 
  typedef std::complex<double> C; 
  
  std::vector<R> MassMatrix {
    10.0,  0.0,
     0.0, 15.0
  };
  
  std::vector<R> StiffnessMatrix {
    1500.0, -500.0, 
    -500.0,  500.0 
  };
  
  R DampingRatio = 0.2; 
  R Omega        = 5.0; 
  
  std::vector<C> expected {
    C ( 1250.0,  3000.0 ), C ( -500.0, -1000.0 ), 
    C ( -500.0, -1000.0 ), C (  125.0,  1000.0 ) 
  };
  
  auto result = 
    SimpleModels::MassSpringDamper::DynamicStiffnessMatrix<R,C> (
      MassMatrix, StiffnessMatrix, 
      DampingRatio, Omega
    );
  
  ASSERT_EQ ( expected.size(), result.size() );

  for ( auto i = 0; i < result.size(); i++ ) {
    EXPECT_EQ ( expected[i].real(), result[i].real() );
    EXPECT_EQ ( expected[i].imag(), result[i].imag() );
  }

} // DynamicStiffnessMatrix, NDOF_2_WithDamping 


TEST ( DynamicStiffnessMatrix, NDOF_3_NoDamping ) {
  
  typedef double               R; 
  typedef std::complex<double> C; 
  
  std::vector<R> MassMatrix {
    20.0,  0.0,  0.0, 
     0.0, 15.0,  0.0, 
     0.0,  0.0, 10.0
  };
  
  std::vector<R> StiffnessMatrix {
    1500.0, -500.0,    0.0, 
    -500.0, 1000.0, -500.0, 
       0.0, -500.0,  500.0 
  };
  
  R DampingRatio = 0.0; 
  R Omega        = 1.0; 
  
  std::vector<C> expected {
    C ( 1480.0, 0.0 ), C ( -500.0, 0.0 ), C (    0.0, 0.0 ), 
    C ( -500.0, 0.0 ), C (  985.0, 0.0 ), C ( -500.0, 0.0 ), 
    C (    0.0, 0.0 ), C ( -500.0, 0.0 ), C (  490.0, 0.0 ) 
  };
  
  auto result = 
    SimpleModels::MassSpringDamper::DynamicStiffnessMatrix<R,C> (
      MassMatrix, StiffnessMatrix, 
      DampingRatio, Omega
    );
  
  ASSERT_EQ ( expected.size(), result.size() );

  for ( auto i = 0; i < result.size(); i++ ) {
    EXPECT_EQ ( expected[i].real(), result[i].real() );
    EXPECT_EQ ( expected[i].imag(), result[i].imag() );
  }

} // DynamicStiffnessMatrix, NDOF_3_NoDamping 


TEST ( DynamicStiffnessMatrix, NDOF_3_WithDamping ) {
  
  typedef double               R; 
  typedef std::complex<double> C; 
  
  std::vector<R> MassMatrix {
    20.0,  0.0,  0.0, 
     0.0, 15.0,  0.0, 
     0.0,  0.0, 10.0
  };
  
  std::vector<R> StiffnessMatrix {
    1500.0, -500.0,    0.0, 
    -500.0, 1000.0, -500.0, 
       0.0, -500.0,  500.0 
  };
  
  R DampingRatio = 0.1; 
  R Omega        = 1.0; 
  
  std::vector<C> expected {
    C ( 1480.0,  300.0 ), C ( -500.0, -100.0 ), C (     0.0,    0.0 ), 
    C ( -500.0, -100.0 ), C (  985.0,  200.0 ), C (  -500.0, -100.0 ), 
    C (    0.0,    0.0 ), C ( -500.0, -100.0 ), C (  490.0,   100.0 ) 
  };
  
  auto result = 
    SimpleModels::MassSpringDamper::DynamicStiffnessMatrix<R,C> (
      MassMatrix, StiffnessMatrix, 
      DampingRatio, Omega
    );
  
  ASSERT_EQ ( expected.size(), result.size() );

  for ( auto i = 0; i < result.size(); i++ ) {
    EXPECT_EQ ( expected[i].real(), result[i].real() );
    EXPECT_EQ ( expected[i].imag(), result[i].imag() );
  }

} // DynamicStiffnessMatrix, NDOF_3_WithDamping 


TEST ( DynamicStiffnessMatrix, NDOF_4_NoDamping ) {
  
  typedef double               R; 
  typedef std::complex<double> C; 
  
  std::vector<R> MassMatrix {
    20.0,  0.0,  0.0,  0.0, 
     0.0, 15.0,  0.0,  0.0, 
     0.0,  0.0, 10.0,  0.0, 
     0.0,  0.0,  0.0, 10.0 
  };
  
  std::vector<R> StiffnessMatrix {
    1500.0, -500.0,    0.0,    0.0, 
    -500.0, 1000.0, -500.0,    0.0, 
       0.0, -500.0,  800.0, -300.0, 
       0.0,    0.0, -300.0,  300.0 
  };
  
  R DampingRatio = 0.0; 
  R Omega        = 1.0; 
  
  std::vector<C> expected {
    C ( 1480.0, 0.0 ), C ( -500.0, 0.0 ), C (    0.0, 0.0 ), C (    0.0, 0.0 ), 
    C ( -500.0, 0.0 ), C (  985.0, 0.0 ), C ( -500.0, 0.0 ), C (    0.0, 0.0 ), 
    C (    0.0, 0.0 ), C ( -500.0, 0.0 ), C (  790.0, 0.0 ), C ( -300.0, 0.0 ),
    C (    0.0, 0.0 ), C (    0.0, 0.0 ), C ( -300.0, 0.0 ), C (  290.0, 0.0 )
  };
  
  auto result = 
    SimpleModels::MassSpringDamper::DynamicStiffnessMatrix<R,C> (
      MassMatrix, StiffnessMatrix, 
      DampingRatio, Omega
    );
  
  ASSERT_EQ ( expected.size(), result.size() );

  for ( auto i = 0; i < result.size(); i++ ) {
    EXPECT_EQ ( expected[i].real(), result[i].real() );
    EXPECT_EQ ( expected[i].imag(), result[i].imag() );
  }

} // DynamicStiffnessMatrix, NDOF_4_NoDamping 


TEST ( DynamicStiffnessMatrix, NDOF_4_WithDamping ) {
  
  typedef double               R; 
  typedef std::complex<double> C; 
  
  std::vector<R> MassMatrix {
    20.0,  0.0,  0.0,  0.0, 
     0.0, 15.0,  0.0,  0.0, 
     0.0,  0.0, 10.0,  0.0, 
     0.0,  0.0,  0.0, 10.0 
  };
  
  std::vector<R> StiffnessMatrix {
    1500.0, -500.0,    0.0,    0.0, 
    -500.0, 1000.0, -500.0,    0.0, 
       0.0, -500.0,  800.0, -300.0, 
       0.0,    0.0, -300.0,  300.0 
  };
  
  R DampingRatio = 0.1; 
  R Omega        = 1.0; 
  
  std::vector<C> expected {

    C ( 1480.0, 300.0 ), C ( -500.0, -100.0 ), 
    C (    0.0,   0.0 ), C (    0.0,    0.0 ), 

    C ( -500.0, -100.0 ), C (  985.0, 200.0 ), 
    C ( -500.0, -100.0 ), C (    0.0,   0.0 ), 

    C (    0.0,   0.0 ), C ( -500.0, -100.0 ), 
    C (  790.0, 160.0 ), C ( -300.0,  -60.0 ),

    C (    0.0,   0.0 ), C (    0.0,  0.0 ), 
    C ( -300.0, -60.0 ), C (  290.0, 60.0 )

  };
  
  auto result = 
    SimpleModels::MassSpringDamper::DynamicStiffnessMatrix<R,C> (
      MassMatrix, StiffnessMatrix, 
      DampingRatio, Omega
    );
  
  ASSERT_EQ ( expected.size(), result.size() );

  for ( auto i = 0; i < result.size(); i++ ) {
    EXPECT_EQ ( expected[i].real(), result[i].real() );
    EXPECT_EQ ( expected[i].imag(), result[i].imag() );
  }

} // DynamicStiffnessMatrix, NDOF_4_WithDamping 


#ifdef DEBUG_SIMPLE_MODELS 


TEST ( DynamicStiffnessMatrix, NegDampingRatio ) {

  typedef double               R; 
  typedef std::complex<double> C; 
  
  std::vector<R> MassMatrix {
    20.0,  0.0,  0.0,  0.0, 
     0.0, 15.0,  0.0,  0.0, 
     0.0,  0.0, 10.0,  0.0, 
     0.0,  0.0,  0.0, 10.0 
  };
  
  std::vector<R> StiffnessMatrix {
    1500.0, -500.0,    0.0,    0.0, 
    -500.0, 1000.0, -500.0,    0.0, 
       0.0, -500.0,  800.0, -300.0, 
       0.0,    0.0, -300.0,  300.0 
  };
  
  R DampingRatio = -0.1; 
  R Omega        = 1.0; 

  bool exception_thrown = false; 

  try {

    auto result = 
      SimpleModels::MassSpringDamper::DynamicStiffnessMatrix<R,C> (
        MassMatrix, StiffnessMatrix, 
        DampingRatio, Omega
      );

  } catch ( const std::exception& e ) {

    EXPECT_STREQ (
      "MassSpringDamper::DynamicStiffnessMatrix: negative DampingRatio", 
      e.what()
    );

    exception_thrown = true;

  }

  EXPECT_TRUE ( exception_thrown ); 

} // DynamicStiffnessMatrix, NegDampingRatio  


TEST ( DynamicStiffnessMatrix, SizeMismatch ) {

  typedef double               R; 
  typedef std::complex<double> C; 
  
  std::vector<R> MassMatrix {
    20.0,  0.0,  0.0,  0.0, 
     0.0, 15.0,  0.0,  0.0, 
     0.0,  0.0, 10.0,  0.0, 
     0.0,  0.0,  0.0, 10.0 
  };
  
  std::vector<R> StiffnessMatrix {
    1500.0, -500.0,    0.0, 
    -500.0, 1000.0, -500.0, 
       0.0, -500.0,  500.0 
  };
  
  R DampingRatio = 0.1; 
  R Omega        = 1.0; 

  bool exception_thrown = false; 

  try {

    auto result = 
      SimpleModels::MassSpringDamper::DynamicStiffnessMatrix<R,C> (
        MassMatrix, StiffnessMatrix, 
        DampingRatio, Omega
      );

  } catch ( const std::exception& e ) {

    EXPECT_STREQ (
      "MassSpringDamper::DynamicStiffnessMatrix: matrices size doesn't match", 
      e.what()
    );

    exception_thrown = true;

  }

  EXPECT_TRUE ( exception_thrown ); 

} // DynamicStiffnessMatrix, SizeMismatch 


#endif // DEBUG_SIMPLE_MODELS 


