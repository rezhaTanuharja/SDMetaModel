/**
  * @file StiffnessMatrix_test.cpp 
  *
  * @brief
  * test functions to compute stiffness matrix for mass-spring-damper models 
  * including errors thrown in DEBUG mode
  *
  * @author
  * Tanuharja R.A. @n 
  * contact: rezha.tanuharja@tum.de 
  */

#include "SimpleModels.hpp" 
#include <gtest/gtest.h> 


TEST ( StiffnessMatrix, NDOF_1 ) {

  std::vector<double> Springs {
    1000.0
  };

  std::vector<double> expected {
    1000.0
  };

  auto result =
    SimpleModels::MassSpringDamper::StiffnessMatrix<double> ( Springs );

  ASSERT_EQ ( expected.size(), result.size() );

  for ( auto i = 0; i < result.size(); i++ ) {
    EXPECT_EQ ( expected[i], result[i] );
  }

} // StiffnessMatrix, NDOF_1 


TEST ( StiffnessMatrix, NDOF_2 ) {

  std::vector<double> Springs {
    1000.0, 500.0 
  };

  std::vector<double> expected {
    1500.0, -500.0, 
    -500.0,  500.0 
  };

  auto result =
    SimpleModels::MassSpringDamper::StiffnessMatrix<double> ( Springs );

  ASSERT_EQ ( expected.size(), result.size() );

  for ( auto i = 0; i < result.size(); i++ ) {
    EXPECT_EQ ( expected[i], result[i] );
  }

} // StiffnessMatrix, NDOF_2 


TEST ( StiffnessMatrix, NDOF_3 ) {

  std::vector<double> Springs {
    1000.0, 500.0, 500.0 
  };

  std::vector<double> expected {
    1500.0, -500.0,    0.0, 
    -500.0, 1000.0, -500.0, 
       0.0, -500.0,  500.0 
  };

  auto result =
    SimpleModels::MassSpringDamper::StiffnessMatrix<double> ( Springs );

  ASSERT_EQ ( expected.size(), result.size() );

  for ( auto i = 0; i < result.size(); i++ ) {
    EXPECT_EQ ( expected[i], result[i] );
  }

} // StiffnessMatrix, NDOF_3 


TEST ( StiffnessMatrix, NDOF_4 ) {

  std::vector<double> Springs {
    1000.0, 500.0, 500.0, 300.0 
  };

  std::vector<double> expected {
    1500.0, -500.0,    0.0,    0.0, 
    -500.0, 1000.0, -500.0,    0.0, 
       0.0, -500.0,  800.0, -300.0, 
       0.0,    0.0, -300.0,  300.0 
  };

  auto result =
    SimpleModels::MassSpringDamper::StiffnessMatrix<double> ( Springs );

  ASSERT_EQ ( expected.size(), result.size() );

  for ( auto i = 0; i < result.size(); i++ ) {
    EXPECT_EQ ( expected[i], result[i] );
  }

} // StiffnessMatrix, NDOF_4 


