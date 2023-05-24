/**
  * @file RandomSampling_test.cpp 
  *
  * @brief
  * test functions to generate random samples from @ref _Normal_ distribution 
  * including errors thrown in DEBUG mode
  *
  * @author
  * Tanuharja R.A. @n 
  * contact: rezha.tanuharja@tum.de 
  */

#include "Statistics.hpp" 
#include <gtest/gtest.h> 


#ifdef DEBUG_STATISTICS 


TEST ( NormalRandomSampling, ZeroPoints ) {
  
  typedef size_t               Z; 
  typedef double               R; 
  typedef std::complex<double> C; 
  
  Z nPoints = 0; 
  Z Dim     = 1; 

  R Mean   =  3.0; 
  R StdDev = 10.0; 

  bool exception_thrown = false; 

  try {
    
    auto result = Statistics::Normal::RandomSampling<Z,R,C> (
      nPoints, Dim,
      Mean, StdDev
    );

  } catch ( const std::exception& e ) {

    EXPECT_STREQ (

      "Normal::RandomSampling: nPoints not in [1,80000] range", 
      e.what() 

    );

    exception_thrown = true; 

  }

  EXPECT_TRUE ( exception_thrown );

} // NormalRandomSampling, ZeroPoints   


TEST ( NormalRandomSampling, NegPoints ) {
  
  typedef size_t               Z; 
  typedef double               R; 
  typedef std::complex<double> C; 
  
  Z nPoints = -1; 
  Z Dim     =  1; 

  R Mean   =  3.0; 
  R StdDev = 10.0; 

  bool exception_thrown = false; 

  try {
    
    auto result = Statistics::Normal::RandomSampling<Z,R,C> (
      nPoints, Dim,
      Mean, StdDev
    );

  } catch ( const std::exception& e ) {

    EXPECT_STREQ (

      "Normal::RandomSampling: nPoints not in [1,80000] range", 
      e.what() 

    );

    exception_thrown = true; 

  }

  EXPECT_TRUE ( exception_thrown );

} // NormalRandomSampling, NegPoints   


TEST ( NormalRandomSampling, HugePoints ) {
  
  typedef size_t               Z; 
  typedef double               R; 
  typedef std::complex<double> C; 
  
  Z nPoints = 200000; 
  Z Dim     =  1; 

  R Mean   =  3.0; 
  R StdDev = 10.0; 

  bool exception_thrown = false; 

  try {
    
    auto result = Statistics::Normal::RandomSampling<Z,R,C> (
      nPoints, Dim,
      Mean, StdDev
    );

  } catch ( const std::exception& e ) {

    EXPECT_STREQ (

      "Normal::RandomSampling: nPoints not in [1,80000] range", 
      e.what() 

    );

    exception_thrown = true; 

  }

  EXPECT_TRUE ( exception_thrown );

} // NormalRandomSampling, HugePoints   


TEST ( NormalRandomSampling, ZeroDim ) {
  
  typedef size_t               Z; 
  typedef double               R; 
  typedef std::complex<double> C; 
  
  Z nPoints = 1; 
  Z Dim     = 0; 

  R Mean   =  3.0; 
  R StdDev = 10.0; 

  bool exception_thrown = false; 

  try {
    
    auto result = Statistics::Normal::RandomSampling<Z,R,C> (
      nPoints, Dim,
      Mean, StdDev
    );

  } catch ( const std::exception& e ) {

    EXPECT_STREQ (

      "Normal::RandomSampling: Dim not in [1,100] range", 
      e.what() 

    );

    exception_thrown = true; 

  }

  EXPECT_TRUE ( exception_thrown );

} // NormalRandomSampling, ZeroDim   


TEST ( NormalRandomSampling, NegDim ) {
  
  typedef size_t               Z; 
  typedef double               R; 
  typedef std::complex<double> C; 
  
  Z nPoints =  1; 
  Z Dim     = -2; 

  R Mean   =  3.0; 
  R StdDev = 10.0; 

  bool exception_thrown = false; 

  try {
    
    auto result = Statistics::Normal::RandomSampling<Z,R,C> (
      nPoints, Dim,
      Mean, StdDev
    );

  } catch ( const std::exception& e ) {

    EXPECT_STREQ (

      "Normal::RandomSampling: Dim not in [1,100] range", 
      e.what() 

    );

    exception_thrown = true; 

  }

  EXPECT_TRUE ( exception_thrown );

} // NormalRandomSampling, NegDim   


TEST ( NormalRandomSampling, HugeDim ) {
  
  typedef size_t               Z; 
  typedef double               R; 
  typedef std::complex<double> C; 
  
  Z nPoints =  1; 
  Z Dim     = 200; 

  R Mean   =  3.0; 
  R StdDev = 10.0; 

  bool exception_thrown = false; 

  try {
    
    auto result = Statistics::Normal::RandomSampling<Z,R,C> (
      nPoints, Dim,
      Mean, StdDev
    );

  } catch ( const std::exception& e ) {

    EXPECT_STREQ (

      "Normal::RandomSampling: Dim not in [1,100] range", 
      e.what() 

    );

    exception_thrown = true; 

  }

  EXPECT_TRUE ( exception_thrown );

} // NormalRandomSampling, HugeDim   


TEST ( NormalRandomSampling, NegStdDev ) {
  
  typedef size_t               Z; 
  typedef double               R; 
  typedef std::complex<double> C; 
  
  Z nPoints =  1; 
  Z Dim     = 2; 

  R Mean   =  3.0; 
  R StdDev = -10.0; 

  bool exception_thrown = false; 

  try {
    
    auto result = Statistics::Normal::RandomSampling<Z,R,C> (
      nPoints, Dim,
      Mean, StdDev
    );

  } catch ( const std::exception& e ) {

    EXPECT_STREQ (

      "Normal::RandomSampling: negative StdDev", 
      e.what() 

    );

    exception_thrown = true; 

  }

  EXPECT_TRUE ( exception_thrown );

} // NormalRandomSampling, HugeDim   


#endif // DEBUG_STATISTICS 

