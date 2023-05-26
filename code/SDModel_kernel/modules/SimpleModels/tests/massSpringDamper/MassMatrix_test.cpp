/**
  * @file MassMatrix_test.cpp 
  *
  * @brief
  * test functions to compute mass matrix for mass-spring-damper models 
  * including errors thrown in DEBUG mode
  *
  * @author
  * Tanuharja R.A. @n 
  * contact: rezha.tanuharja@tum.de 
  */

#include "SimpleModels.hpp" 
#include <gtest/gtest.h> 


TEST ( MassMatrix, NDOF_1 ) {

  std::vector<double> Masses {
    15.0
  };

  std::vector<double> expected {
     15.0
  };

  auto result = 
    SimpleModels::MassSpringDamper::MassMatrix<double> ( Masses );

  ASSERT_EQ ( expected.size(), result.size() );

  for ( auto i = 0; i < result.size(); i++ ) {
    EXPECT_EQ ( expected[i], result[i] );
  }

} // MassMatrix, NDOF_1 


TEST ( MassMatrix, NDOF_2 ) {

  std::vector<double> Masses {
    10.0, 15.0
  };

  std::vector<double> expected {
    10.0,  0.0,
     0.0, 15.0
  };

  auto result = 
    SimpleModels::MassSpringDamper::MassMatrix<double> ( Masses );

  ASSERT_EQ ( expected.size(), result.size() );

  for ( auto i = 0; i < result.size(); i++ ) {
    EXPECT_EQ ( expected[i], result[i] );
  }

} // MassMatrix, NDOF_2 


TEST ( MassMatrix, NDOF_3 ) {

  std::vector<double> Masses {
    20.0, 15.0, 10.0 
  };

  std::vector<double> expected {
    20.0,  0.0,  0.0, 
     0.0, 15.0,  0.0, 
     0.0,  0.0, 10.0
  };

  auto result = 
    SimpleModels::MassSpringDamper::MassMatrix<double> ( Masses );

  ASSERT_EQ ( expected.size(), result.size() );

  for ( auto i = 0; i < result.size(); i++ ) {
    EXPECT_EQ ( expected[i], result[i] );
  }

} // MassMatrix, NDOF_3 


TEST ( MassMatrix, NDOF_4 ) {

  std::vector<double> Masses {
    20.0, 15.0, 10.0, 10.0 
  };

  std::vector<double> expected {
    20.0,  0.0,  0.0,  0.0, 
     0.0, 15.0,  0.0,  0.0, 
     0.0,  0.0, 10.0,  0.0, 
     0.0,  0.0,  0.0, 10.0 
  };

  auto result = 
    SimpleModels::MassSpringDamper::MassMatrix<double> ( Masses );

  ASSERT_EQ ( expected.size(), result.size() );

  for ( auto i = 0; i < result.size(); i++ ) {
    EXPECT_EQ ( expected[i], result[i] );
  }

} // MassMatrix, NDOF_4 


