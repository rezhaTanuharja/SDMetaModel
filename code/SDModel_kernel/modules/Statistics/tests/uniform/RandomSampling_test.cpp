/**
  * @file RandomSampling_test.cpp 
  *
  * @brief
  * test functions to generate random samples from @ref _Uniform_ distribution 
  * including errors thrown in DEBUG mode
  *
  * @author
  * Tanuharja R.A. @n 
  * contact: rezha.tanuharja@tum.de 
  */

#include "Statistics.hpp" 
#include <gtest/gtest.h> 


TEST ( UniformRandomSampling, StdUniform_MultiPtsMultiDims ) {
  
  typedef size_t               Z; 
  typedef double               R; 
  typedef std::complex<double> C; 
  
  Z nPoints = 40; 
  Z Dim     = 3; 

  R LowerBound = 0.0; 
  R UpperBound = 1.0; 

  auto result = Statistics::Uniform::RandomSampling<Z,R,C> (
    nPoints, Dim,
    LowerBound, UpperBound
  );

  ASSERT_EQ ( result.size(), nPoints * Dim );

  std::vector<bool> CheckSample ( nPoints * Dim );

  for ( auto i = 0; i < result.size(); i++ ) {
    CheckSample[i] =
       result[i].real() >= LowerBound && 
       result[i].real() <= UpperBound;
  }

  EXPECT_TRUE (
    std::all_of (
      CheckSample.begin(), 
      CheckSample.end(),

      []( const auto m ) { return m; }
    )
  );

} // UniformRandomSampling, StdUniform_MultiPtsMultiDims 


TEST ( UniformRandomSampling, StdUniform_SinglePtsMultiDims ) {
  
  typedef size_t               Z; 
  typedef double               R; 
  typedef std::complex<double> C; 
  
  Z nPoints = 1; 
  Z Dim     = 3; 

  R LowerBound = 0.0; 
  R UpperBound = 1.0; 

  auto result = Statistics::Uniform::RandomSampling<Z,R,C> (
    nPoints, Dim,
    LowerBound, UpperBound
  );

  ASSERT_EQ ( result.size(), nPoints * Dim );

  std::vector<bool> CheckSample ( nPoints * Dim );

  for ( auto i = 0; i < result.size(); i++ ) {
    CheckSample[i] =
       result[i].real() >= LowerBound && 
       result[i].real() <= UpperBound;
  }

  EXPECT_TRUE (
    std::all_of (
      CheckSample.begin(), 
      CheckSample.end(),

      []( const auto m ) { return m; }
    )
  );

} // UniformRandomSampling, StdUniform_SinglePtsMultiDims 


TEST ( UniformRandomSampling, StdUniform_MultiPtsSingleDims ) {
  
  typedef size_t               Z; 
  typedef double               R; 
  typedef std::complex<double> C; 
  
  Z nPoints = 100; 
  Z Dim     = 1; 

  R LowerBound = 0.0; 
  R UpperBound = 1.0; 

  auto result = Statistics::Uniform::RandomSampling<Z,R,C> (
    nPoints, Dim,
    LowerBound, UpperBound
  );

  ASSERT_EQ ( result.size(), nPoints * Dim );

  std::vector<bool> CheckSample ( nPoints * Dim );

  for ( auto i = 0; i < result.size(); i++ ) {
    CheckSample[i] =
       result[i].real() >= LowerBound && 
       result[i].real() <= UpperBound;
  }

  EXPECT_TRUE (
    std::all_of (
      CheckSample.begin(), 
      CheckSample.end(),

      []( const auto m ) { return m; }
    )
  );

} // UniformRandomSampling, StdUniform_MultiPtsSingleDims 


TEST ( UniformRandomSampling, StdUniform_SinglePtsSingleDims ) {
  
  typedef size_t               Z; 
  typedef double               R; 
  typedef std::complex<double> C; 
  
  Z nPoints = 1; 
  Z Dim     = 1; 

  R LowerBound = 0.0; 
  R UpperBound = 1.0; 

  auto result = Statistics::Uniform::RandomSampling<Z,R,C> (
    nPoints, Dim,
    LowerBound, UpperBound
  );

  ASSERT_EQ ( result.size(), nPoints * Dim );

  std::vector<bool> CheckSample ( nPoints * Dim );

  for ( auto i = 0; i < result.size(); i++ ) {
    CheckSample[i] =
       result[i].real() >= LowerBound && 
       result[i].real() <= UpperBound;
  }

  EXPECT_TRUE (
    std::all_of (
      CheckSample.begin(), 
      CheckSample.end(),

      []( const auto m ) { return m; }
    )
  );

} // UniformRandomSampling, StdUniform_SinglePtsSingleDims 


TEST ( UniformRandomSampling, MultiPtsMultiDims ) {
  
  typedef size_t               Z; 
  typedef double               R; 
  typedef std::complex<double> C; 
  
  Z nPoints = 40; 
  Z Dim     = 3; 

  R LowerBound = -1.0; 
  R UpperBound =  7.0; 

  auto result = Statistics::Uniform::RandomSampling<Z,R,C> (
    nPoints, Dim,
    LowerBound, UpperBound
  );

  ASSERT_EQ ( result.size(), nPoints * Dim );

  std::vector<bool> CheckSample ( nPoints * Dim );

  for ( auto i = 0; i < result.size(); i++ ) {
    CheckSample[i] =
       result[i].real() >= LowerBound && 
       result[i].real() <= UpperBound;
  }

  EXPECT_TRUE (
    std::all_of (
      CheckSample.begin(), 
      CheckSample.end(),

      []( const auto m ) { return m; }
    )
  );

} // UniformRandomSampling, MultiPtsMultiDims  


TEST ( UniformRandomSampling, SinglePtsMultiDims ) {
  
  typedef size_t               Z; 
  typedef double               R; 
  typedef std::complex<double> C; 
  
  Z nPoints = 1; 
  Z Dim     = 3; 

  R LowerBound = -1.0; 
  R UpperBound =  7.0; 

  auto result = Statistics::Uniform::RandomSampling<Z,R,C> (
    nPoints, Dim,
    LowerBound, UpperBound
  );

  ASSERT_EQ ( result.size(), nPoints * Dim );

  std::vector<bool> CheckSample ( nPoints * Dim );

  for ( auto i = 0; i < result.size(); i++ ) {
    CheckSample[i] =
       result[i].real() >= LowerBound && 
       result[i].real() <= UpperBound;
  }

  EXPECT_TRUE (
    std::all_of (
      CheckSample.begin(), 
      CheckSample.end(),

      []( const auto m ) { return m; }
    )
  );

} // UniformRandomSampling, SinglePtsMultiDims  


TEST ( UniformRandomSampling, MultiPtsSingleDims ) {
  
  typedef size_t               Z; 
  typedef double               R; 
  typedef std::complex<double> C; 
  
  Z nPoints = 100; 
  Z Dim     = 1; 

  R LowerBound = 3.0; 
  R UpperBound = 7.0; 

  auto result = Statistics::Uniform::RandomSampling<Z,R,C> (
    nPoints, Dim,
    LowerBound, UpperBound
  );

  ASSERT_EQ ( result.size(), nPoints * Dim );

  std::vector<bool> CheckSample ( nPoints * Dim );

  for ( auto i = 0; i < result.size(); i++ ) {
    CheckSample[i] =
       result[i].real() >= LowerBound && 
       result[i].real() <= UpperBound;
  }

  EXPECT_TRUE (
    std::all_of (
      CheckSample.begin(), 
      CheckSample.end(),

      []( const auto m ) { return m; }
    )
  );

} // UniformRandomSampling, MultiPtsSingleDims   


TEST ( UniformRandomSampling, SinglePtsSingleDims ) {
  
  typedef size_t               Z; 
  typedef double               R; 
  typedef std::complex<double> C; 
  
  Z nPoints = 1; 
  Z Dim     = 1; 

  R LowerBound =  3.0; 
  R UpperBound = 10.0; 

  auto result = Statistics::Uniform::RandomSampling<Z,R,C> (
    nPoints, Dim,
    LowerBound, UpperBound
  );

  ASSERT_EQ ( result.size(), nPoints * Dim );

  std::vector<bool> CheckSample ( nPoints * Dim );

  for ( auto i = 0; i < result.size(); i++ ) {
    CheckSample[i] =
       result[i].real() >= LowerBound && 
       result[i].real() <= UpperBound;
  }

  EXPECT_TRUE (
    std::all_of (
      CheckSample.begin(), 
      CheckSample.end(),

      []( const auto m ) { return m; }
    )
  );

} // UniformRandomSampling, SinglePtsSingleDims   


#ifdef DEBUG_STATISTICS 


TEST ( UniformRandomSampling, ZeroPoints ) {
  
  typedef size_t               Z; 
  typedef double               R; 
  typedef std::complex<double> C; 
  
  Z nPoints = 0; 
  Z Dim     = 1; 

  R LowerBound =  3.0; 
  R UpperBound = 10.0; 

  bool exception_thrown = false; 

  try {
    
    auto result = Statistics::Uniform::RandomSampling<Z,R,C> (
      nPoints, Dim,
      LowerBound, UpperBound
    );

  } catch ( const std::exception& e ) {

    EXPECT_STREQ (

      "Uniform::RandomSampling: nPoints not in [1,80000] range", 
      e.what() 

    );

    exception_thrown = true; 

  }

  EXPECT_TRUE ( exception_thrown );

} // UniformRandomSampling, ZeroPoints   


TEST ( UniformRandomSampling, NegPoints ) {
  
  typedef size_t               Z; 
  typedef double               R; 
  typedef std::complex<double> C; 
  
  Z nPoints = -2; 
  Z Dim     = 1; 

  R LowerBound =  3.0; 
  R UpperBound = 10.0; 

  bool exception_thrown = false; 

  try {
    
    auto result = Statistics::Uniform::RandomSampling<Z,R,C> (
      nPoints, Dim,
      LowerBound, UpperBound
    );

  } catch ( const std::exception& e ) {

    EXPECT_STREQ (

      "Uniform::RandomSampling: nPoints not in [1,80000] range", 
      e.what() 

    );

    exception_thrown = true; 

  }

  EXPECT_TRUE ( exception_thrown );

} // UniformRandomSampling, NegPoints 


TEST ( UniformRandomSampling, HugePoints ) {
  
  typedef size_t               Z; 
  typedef double               R; 
  typedef std::complex<double> C; 
  
  Z nPoints = 200000; 
  Z Dim     = 1; 

  R LowerBound =  3.0; 
  R UpperBound = 10.0; 

  bool exception_thrown = false; 

  try {
    
    auto result = Statistics::Uniform::RandomSampling<Z,R,C> (
      nPoints, Dim,
      LowerBound, UpperBound
    );

  } catch ( const std::exception& e ) {

    EXPECT_STREQ (

      "Uniform::RandomSampling: nPoints not in [1,80000] range", 
      e.what() 

    );

    exception_thrown = true; 

  }

  EXPECT_TRUE ( exception_thrown );

} // UniformRandomSampling, HugePoints 


TEST ( UniformRandomSampling, ZeroDim ) {
  
  typedef size_t               Z; 
  typedef double               R; 
  typedef std::complex<double> C; 
  
  Z nPoints = 1; 
  Z Dim     = 0; 

  R LowerBound =  3.0; 
  R UpperBound = 10.0; 

  bool exception_thrown = false; 

  try {
    
    auto result = Statistics::Uniform::RandomSampling<Z,R,C> (
      nPoints, Dim,
      LowerBound, UpperBound
    );

  } catch ( const std::exception& e ) {

    EXPECT_STREQ (

      "Uniform::RandomSampling: Dim not in [1,100] range", 
      e.what() 

    );

    exception_thrown = true; 

  }

  EXPECT_TRUE ( exception_thrown );

} // UniformRandomSampling, ZeroDim    


TEST ( UniformRandomSampling, NegDim ) {
  
  typedef size_t               Z; 
  typedef double               R; 
  typedef std::complex<double> C; 
  
  Z nPoints =  1; 
  Z Dim     = -5; 

  R LowerBound =  3.0; 
  R UpperBound = 10.0; 

  bool exception_thrown = false; 

  try {
    
    auto result = Statistics::Uniform::RandomSampling<Z,R,C> (
      nPoints, Dim,
      LowerBound, UpperBound
    );

  } catch ( const std::exception& e ) {

    EXPECT_STREQ (

      "Uniform::RandomSampling: Dim not in [1,100] range", 
      e.what() 

    );

    exception_thrown = true; 

  }

  EXPECT_TRUE ( exception_thrown );

} // UniformRandomSampling, NegDim    


TEST ( UniformRandomSampling, HugeDim ) {
  
  typedef size_t               Z; 
  typedef double               R; 
  typedef std::complex<double> C; 
  
  Z nPoints =  1; 
  Z Dim     = 500; 

  R LowerBound =  3.0; 
  R UpperBound = 10.0; 

  bool exception_thrown = false; 

  try {
    
    auto result = Statistics::Uniform::RandomSampling<Z,R,C> (
      nPoints, Dim,
      LowerBound, UpperBound
    );

  } catch ( const std::exception& e ) {

    EXPECT_STREQ (

      "Uniform::RandomSampling: Dim not in [1,100] range", 
      e.what() 

    );

    exception_thrown = true; 

  }

  EXPECT_TRUE ( exception_thrown );

} // UniformRandomSampling, HugeDim    


TEST ( UniformRandomSampling, FlippedBounds ) {
  
  typedef size_t               Z; 
  typedef double               R; 
  typedef std::complex<double> C; 
  
  Z nPoints = 3; 
  Z Dim     = 5; 

  R LowerBound = 30.0; 
  R UpperBound = 10.0; 

  bool exception_thrown = false; 

  try {
    
    auto result = Statistics::Uniform::RandomSampling<Z,R,C> (
      nPoints, Dim,
      LowerBound, UpperBound
    );

  } catch ( const std::exception& e ) {

    EXPECT_STREQ (

      "Uniform::RandomSampling: LowerBound >= UpperBound", 
      e.what() 

    );

    exception_thrown = true; 

  }

  EXPECT_TRUE ( exception_thrown );

} // UniformRandomSampling, FlippedBounds 


#endif // DEBUG_STATISTICS 

