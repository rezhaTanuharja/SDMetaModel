/**
  * @file MassSpringDamper_DMC_test.cpp 
  *
  * @brief
  * test support for direct MC simulation of mass-spring-damper model 
  * including errors thrown in DEBUG mode
  *
  * @author
  * Tanuharja R.A. @n 
  * contact: rezha.tanuharja@tum.de 
  */

#include "NIRPCE.hpp" 
#include <gtest/gtest.h> 


TEST ( NIRPCE_Train, RealNDOF1_DIM1 ) {

  typedef Stochastic::Z Z; 
  typedef Stochastic::R R; 
  typedef Stochastic::C C; 

  Z nDOFs = 1; 
  Z Dim   = 1; 

  std::vector<Z> NumIndices {
    0
  };

  std::vector<Z> DenIndices {
    0, 1 
  };

  Stochastic::Surrogate::NIRPCE Model ( nDOFs, Dim );

  Model.SetNumIndices ( NumIndices ); 
  Model.SetDenIndices ( DenIndices ); 

  std::vector<C> InputVars {
    C ( 0.1, 0.0 ), 
    C ( 0.4, 0.0 ), 
    C ( 0.3, 0.0 ) 
  };

  std::vector<C> Responses ( InputVars.size() );

  std::transform (
    InputVars.begin(), 
    InputVars.end(), 

    Responses.begin(), 

    [=]( const auto m ) {
      return C ( 1.0,0.0 ) / m;
    }
  );

  ASSERT_NO_THROW (
    Model.Train ( InputVars, Responses );
  );


  #ifdef DEBUG_RSMSD 

  auto numCoeffs = Model.NumCoeffs(); 
  auto denCoeffs = Model.DenCoeffs(); 

  R tolerance = 1e-4; 

  for ( auto i = 0; i < numCoeffs.size(); i++ ) {
    EXPECT_NEAR ( numCoeffs[i].imag(), 0.0, tolerance );
  }

  for ( auto i = 0; i < denCoeffs.size(); i++ ) {
    EXPECT_NEAR ( denCoeffs[i].imag(), 0.0, tolerance );
  }

  R RelativeError = std::abs ( 
    (numCoeffs[0].real() - denCoeffs[1].real()) / numCoeffs[0].real()
  );

  EXPECT_NEAR ( RelativeError, 0.0 , tolerance );

  #endif // DEBUG_RSMSD 


} // NIRPCE_Train, RealNDOF1_DIM1  


TEST ( NIRPCE_Train, ComplexNDOF1_DIM1 ) {

  typedef Stochastic::Z Z; 
  typedef Stochastic::R R; 
  typedef Stochastic::C C; 

  Z nDOFs = 1; 
  Z Dim   = 1; 

  std::vector<Z> NumIndices {
    0
  };

  std::vector<Z> DenIndices {
    0, 1 
  };

  Stochastic::Surrogate::NIRPCE Model ( nDOFs, Dim );

  Model.SetNumIndices ( NumIndices ); 
  Model.SetDenIndices ( DenIndices ); 

  std::vector<C> InputVars {
    C ( 0.1,  2.0 ), 
    C ( 0.4, -1.0 ), 
    C ( 0.3,  0.0 ) 
  };

  std::vector<C> Responses ( InputVars.size() );

  std::transform (
    InputVars.begin(), 
    InputVars.end(), 

    Responses.begin(), 

    [=]( const auto m ) {
      return C ( 1.0,0.0 ) / m;
    }
  );

  ASSERT_NO_THROW (
    Model.Train ( InputVars, Responses );
  );


  #ifdef DEBUG_RSMSD 

  auto numCoeffs = Model.NumCoeffs(); 
  auto denCoeffs = Model.DenCoeffs(); 

  R tolerance = 1e-4; 

  EXPECT_NEAR ( denCoeffs[0].real(), 0.0, tolerance );
  EXPECT_NEAR ( denCoeffs[0].imag(), 0.0, tolerance );

  EXPECT_NEAR ( numCoeffs[0].real(), denCoeffs[1].real(), tolerance );
  EXPECT_NEAR ( numCoeffs[0].imag(), denCoeffs[1].imag(), tolerance );

  #endif // DEBUG_RSMSD 


} // NIRPCE_Train, ComplexNDOF1_DIM1 


TEST ( NIRPCE_ComputeResponse, RealNDOF1_DIM1 ) {

  typedef Stochastic::Z Z; 
  typedef Stochastic::R R; 
  typedef Stochastic::C C; 

  Z nDOFs = 1; 
  Z Dim   = 1; 

  std::vector<Z> NumIndices {
    0
  };

  std::vector<Z> DenIndices {
    0, 1 
  };

  Stochastic::Surrogate::NIRPCE Model ( nDOFs, Dim );

  Model.SetNumIndices ( NumIndices ); 
  Model.SetDenIndices ( DenIndices ); 

  std::vector<C> InputVars {
    C ( 0.1, 0.0 ), 
    C ( 0.4, 0.0 ), 
    C ( 0.3, 0.0 ) 
  };

  std::vector<C> Responses ( InputVars.size() );

  std::transform (
    InputVars.begin(), 
    InputVars.end(), 

    Responses.begin(), 

    [=]( const auto m ) {
      return C ( 1.0,0.0 ) / m;
    }
  );

  ASSERT_NO_THROW (
    Model.Train ( InputVars, Responses );
  );

  auto result = Model.ComputeResponse ( InputVars ); 

  R tolerance = 1e-5; 

  ASSERT_EQ ( result.size(), Responses.size() );

  for ( auto i = 0; i < result.size(); i++ ) {
    EXPECT_NEAR ( result[i].real(), Responses[i].real(), tolerance );
    EXPECT_NEAR ( result[i].imag(), Responses[i].imag(), tolerance );
  }

} // NIRPCE_ComputeResponse, RealNDOF1_DIM1  


TEST ( NIRPCE_ComputeResponse, ComplexNDOF1_DIM1 ) {

  typedef Stochastic::Z Z; 
  typedef Stochastic::R R; 
  typedef Stochastic::C C; 

  Z nDOFs = 1; 
  Z Dim   = 1; 

  std::vector<Z> NumIndices {
    0
  };

  std::vector<Z> DenIndices {
    0, 1 
  };

  Stochastic::Surrogate::NIRPCE Model ( nDOFs, Dim );

  Model.SetNumIndices ( NumIndices ); 
  Model.SetDenIndices ( DenIndices ); 

  std::vector<C> InputVars {
    C ( 0.1, 2.0 ), 
    C ( 0.4, 0.5 ), 
    C ( 0.0, 0.5 ), 
    C ( 0.4, 0.0 ), 
    C (-0.4, 0.5 ), 
    C ( 1.4, 1.5 ), 
    C ( 0.8, 2.5 ), 
    C ( 0.3,-1.0 ) 
  };

  std::vector<C> Responses ( InputVars.size() );

  std::transform (
    InputVars.begin(), 
    InputVars.end(), 

    Responses.begin(), 

    [=]( const auto m ) {
      return C ( 1.0,0.0 ) / m;
    }
  );

  ASSERT_NO_THROW (
    Model.Train ( InputVars, Responses );
  );

  auto result = Model.ComputeResponse ( InputVars ); 

  R tolerance = 1e-6; 

  ASSERT_EQ ( result.size(), Responses.size() );

  for ( auto i = 0; i < result.size(); i++ ) {
    EXPECT_NEAR ( result[i].real(), Responses[i].real(), tolerance );
    EXPECT_NEAR ( result[i].imag(), Responses[i].imag(), tolerance );
  }

} // NIRPCE_ComputeResponse, ComplexNDOF1_DIM1  


TEST ( NIRPCE_Train, RealNDOF1_DIM2 ) {

  typedef Stochastic::Z Z; 
  typedef Stochastic::R R; 
  typedef Stochastic::C C; 

  Z nDOFs = 1; 
  Z Dim   = 2; 

  std::vector<Z> NumIndices {
    0, 0 
  };

  std::vector<Z> DenIndices {
    0, 0,
    1, 0, 
    0, 1 
  };

  Stochastic::Surrogate::NIRPCE Model ( nDOFs, Dim );

  Model.SetNumIndices ( NumIndices ); 
  Model.SetDenIndices ( DenIndices ); 

  std::vector<C> InputVars {
    C ( 0.1, 0.0 ), 
    C ( 0.4, 0.0 ), 
    C ( 1.0, 0.0 ), 
    C ( 2.4, 0.0 ), 
    C (-0.4, 0.0 ), 
    C ( 1.1, 0.0 ), 
    C ( 0.8, 0.0 ), 
    C ( 0.3, 0.0 ) 
  };

  std::vector<C> Responses ( InputVars.size() / Dim );

  for ( auto i = 0; i < Responses.size(); i++ ) {

    Responses[i] = C ( 1.0, 0.0 ) / (
      InputVars[2*i  ] +
      InputVars[2*i+1]
    );

  }

  ASSERT_NO_THROW (
    Model.Train ( InputVars, Responses );
  );


  #ifdef DEBUG_RSMSD 

  auto numCoeffs = Model.NumCoeffs(); 
  auto denCoeffs = Model.DenCoeffs(); 

  R tolerance = 1e-4; 

  EXPECT_NEAR ( numCoeffs[0].real(), denCoeffs[1].real(), tolerance );
  EXPECT_NEAR ( numCoeffs[0].real(), denCoeffs[2].real(), tolerance );
  EXPECT_NEAR ( denCoeffs[1].real(), denCoeffs[2].real(), tolerance );

  EXPECT_NEAR ( numCoeffs[0].imag(), denCoeffs[1].imag(), tolerance );
  EXPECT_NEAR ( numCoeffs[0].imag(), denCoeffs[2].imag(), tolerance );
  EXPECT_NEAR ( denCoeffs[1].imag(), denCoeffs[2].imag(), tolerance );

  EXPECT_NEAR ( denCoeffs[0].real(), 0.0, tolerance );
  EXPECT_NEAR ( denCoeffs[0].imag(), 0.0, tolerance );

  #endif // DEBUG_RSMSD 


} // NIRPCE_Train, RealNDOF1_DIM2  


TEST ( NIRPCE_Train, ComplexNDOF1_DIM2 ) {

  typedef Stochastic::Z Z; 
  typedef Stochastic::R R; 
  typedef Stochastic::C C; 

  Z nDOFs = 1; 
  Z Dim   = 2; 

  std::vector<Z> NumIndices {
    0, 0 
  };

  std::vector<Z> DenIndices {
    0, 0,
    1, 0, 
    0, 1 
  };

  Stochastic::Surrogate::NIRPCE Model ( nDOFs, Dim );

  Model.SetNumIndices ( NumIndices ); 
  Model.SetDenIndices ( DenIndices ); 

  std::vector<C> InputVars {
    C ( 0.1, 1.0 ), 
    C ( 0.4, 0.6 ), 
    C ( 1.0, 3.0 ), 
    C ( 2.4, 2.0 ), 
    C (-0.4, 2.0 ), 
    C ( 1.1, 0.2 ), 
    C ( 0.8, 0.2 ), 
    C ( 0.3,-0.1 ) 
  };

  std::vector<C> Responses ( InputVars.size() / Dim );

  for ( auto i = 0; i < Responses.size(); i++ ) {

    Responses[i] = C ( 1.0, 0.0 ) / (
      InputVars[2*i  ] +
      InputVars[2*i+1]
    );

  }

  ASSERT_NO_THROW (
    Model.Train ( InputVars, Responses );
  );


  #ifdef DEBUG_RSMSD 

  auto numCoeffs = Model.NumCoeffs(); 
  auto denCoeffs = Model.DenCoeffs(); 

  R tolerance = 1e-4; 

  EXPECT_NEAR ( numCoeffs[0].real(), denCoeffs[1].real(), tolerance );
  EXPECT_NEAR ( numCoeffs[0].real(), denCoeffs[2].real(), tolerance );
  EXPECT_NEAR ( denCoeffs[1].real(), denCoeffs[2].real(), tolerance );

  EXPECT_NEAR ( numCoeffs[0].imag(), denCoeffs[1].imag(), tolerance );
  EXPECT_NEAR ( numCoeffs[0].imag(), denCoeffs[2].imag(), tolerance );
  EXPECT_NEAR ( denCoeffs[1].imag(), denCoeffs[2].imag(), tolerance );

  EXPECT_NEAR ( denCoeffs[0].real(), 0.0, tolerance );
  EXPECT_NEAR ( denCoeffs[0].imag(), 0.0, tolerance );

  #endif // DEBUG_RSMSD 


} // NIRPCE_Train, ComplexNDOF1_DIM2  


TEST ( NIRPCE_ComputeResponse, RealNDOF1_DIM2 ) {

  typedef Stochastic::Z Z; 
  typedef Stochastic::R R; 
  typedef Stochastic::C C; 

  Z nDOFs = 1; 
  Z Dim   = 2; 

  std::vector<Z> NumIndices {
    0, 0 
  };

  std::vector<Z> DenIndices {
    0, 0,
    1, 0, 
    0, 1 
  };

  Stochastic::Surrogate::NIRPCE Model ( nDOFs, Dim );

  Model.SetNumIndices ( NumIndices ); 
  Model.SetDenIndices ( DenIndices ); 

  std::vector<C> InputVars {
    C ( 0.1, 0.0 ), 
    C ( 0.4, 0.0 ), 
    C ( 1.0, 0.0 ), 
    C ( 2.4, 0.0 ), 
    C (-0.4, 0.0 ), 
    C ( 1.1, 0.0 ), 
    C ( 0.8, 0.0 ), 
    C ( 0.3, 0.0 ) 
  };

  std::vector<C> Responses ( InputVars.size() / Dim );

  for ( auto i = 0; i < Responses.size(); i++ ) {

    Responses[i] = C ( 1.0, 0.0 ) / (
      InputVars[2*i  ] +
      InputVars[2*i+1]
    );

  }

  ASSERT_NO_THROW (
    Model.Train ( InputVars, Responses );
  );

  auto result = Model.ComputeResponse ( InputVars );

  R tolerance = 1e-5; 

  ASSERT_EQ ( result.size(), Responses.size() );

  for ( auto i = 0; i < result.size(); i++ ) {
    EXPECT_NEAR ( result[i].real(), Responses[i].real(), tolerance );
    EXPECT_NEAR ( result[i].imag(), Responses[i].imag(), tolerance );
  }

} // NIRPCE_ComputeResponse, RealNDOF1_DIM2  


TEST ( NIRPCE_ComputeResponse, ComplexNDOF1_DIM2 ) {

  typedef Stochastic::Z Z; 
  typedef Stochastic::R R; 
  typedef Stochastic::C C; 

  Z nDOFs = 1; 
  Z Dim   = 2; 

  std::vector<Z> NumIndices {
    0, 0 
  };

  std::vector<Z> DenIndices {
    0, 0,
    1, 0, 
    0, 1 
  };

  Stochastic::Surrogate::NIRPCE Model ( nDOFs, Dim );

  Model.SetNumIndices ( NumIndices ); 
  Model.SetDenIndices ( DenIndices ); 

  std::vector<C> InputVars {
    C ( 0.1, 1.0 ), 
    C ( 0.4, 0.6 ), 
    C ( 1.0, 3.0 ), 
    C ( 2.4, 2.0 ), 
    C (-0.4, 2.0 ), 
    C ( 1.1, 0.2 ), 
    C ( 0.8, 0.2 ), 
    C ( 0.3,-1.0 ) 
  };

  std::vector<C> Responses ( InputVars.size() / Dim );

  for ( auto i = 0; i < Responses.size(); i++ ) {

    Responses[i] = C ( 1.0, 0.0 ) / (
      InputVars[2*i  ] +
      InputVars[2*i+1]
    );

  }

  ASSERT_NO_THROW (
    Model.Train ( InputVars, Responses );
  );

  auto result = Model.ComputeResponse ( InputVars );

  R tolerance = 1e-5; 

  ASSERT_EQ ( result.size(), Responses.size() );

  for ( auto i = 0; i < result.size(); i++ ) {
    EXPECT_NEAR ( result[i].real(), Responses[i].real(), tolerance );
    EXPECT_NEAR ( result[i].imag(), Responses[i].imag(), tolerance );
  }

} // NIRPCE_ComputeResponse, ComplexNDOF1_DIM2  


TEST ( NIRPCE_Train, RealNDOF2_DIM2 ) {

  typedef Stochastic::Z Z; 
  typedef Stochastic::R R; 
  typedef Stochastic::C C; 

  Z nDOFs = 2; 
  Z Dim   = 2; 

  std::vector<Z> NumIndices {
    0, 0 
  };

  std::vector<Z> DenIndices {
    0, 0,
    1, 0, 
    0, 1 
  };

  Stochastic::Surrogate::NIRPCE Model ( nDOFs, Dim );

  Model.SetNumIndices ( NumIndices ); 
  Model.SetDenIndices ( DenIndices ); 

  std::vector<C> InputVars {
    C ( 0.1, 0.0 ), 
    C ( 0.4, 0.0 ), 
    C ( 1.0, 0.0 ), 
    C ( 2.4, 0.0 ), 
    C (-0.4, 0.0 ), 
    C ( 1.1, 0.0 ), 
    C ( 0.8, 0.0 ), 
    C ( 0.3, 0.0 ) 
  };

  std::vector<C> Responses ( nDOFs * InputVars.size() / Dim );

  for ( auto i = 0; i < Responses.size(); i++ ) {

    Responses[2*i  ] = C ( 1.0, 0.0 ) / (
      InputVars[2*i  ] +
      InputVars[2*i+1]
    );

    Responses[2*i+1] = C ( 1.0, 0.0 ) / (
      InputVars[2*i  ]
    );

  }

  ASSERT_NO_THROW (
    Model.Train ( InputVars, Responses );
  );


  #ifdef DEBUG_RSMSD 

  auto numCoeffs = Model.NumCoeffs(); 
  auto denCoeffs = Model.DenCoeffs(); 

  R tolerance = 1e-4; 

  EXPECT_NEAR ( numCoeffs[0].real(), denCoeffs[1].real(), tolerance );
  EXPECT_NEAR ( numCoeffs[0].real(), denCoeffs[2].real(), tolerance );
  EXPECT_NEAR ( denCoeffs[1].real(), denCoeffs[2].real(), tolerance );

  EXPECT_NEAR ( numCoeffs[0].imag(), denCoeffs[1].imag(), tolerance );
  EXPECT_NEAR ( numCoeffs[0].imag(), denCoeffs[2].imag(), tolerance );
  EXPECT_NEAR ( denCoeffs[1].imag(), denCoeffs[2].imag(), tolerance );

  EXPECT_NEAR ( denCoeffs[0].real(), 0.0, tolerance );
  EXPECT_NEAR ( denCoeffs[0].imag(), 0.0, tolerance );

  EXPECT_NEAR ( numCoeffs[1].real(), denCoeffs[4].real(), tolerance );
  EXPECT_NEAR ( numCoeffs[1].imag(), denCoeffs[4].imag(), tolerance );

  EXPECT_NEAR ( denCoeffs[3].real(), 0.0, tolerance );
  EXPECT_NEAR ( denCoeffs[3].imag(), 0.0, tolerance );

  EXPECT_NEAR ( denCoeffs[5].real(), 0.0, tolerance );
  EXPECT_NEAR ( denCoeffs[5].imag(), 0.0, tolerance );

  #endif // DEBUG_RSMSD 


} // NIRPCE_Train, RealNDOF2_DIM2  


TEST ( NIRPCE_Train, ComplexNDOF2_DIM2 ) {

  typedef Stochastic::Z Z; 
  typedef Stochastic::R R; 
  typedef Stochastic::C C; 

  Z nDOFs = 2; 
  Z Dim   = 2; 

  std::vector<Z> NumIndices {
    0, 0 
  };

  std::vector<Z> DenIndices {
    0, 0,
    1, 0, 
    0, 1 
  };

  Stochastic::Surrogate::NIRPCE Model ( nDOFs, Dim );

  Model.SetNumIndices ( NumIndices ); 
  Model.SetDenIndices ( DenIndices ); 

  std::vector<C> InputVars {
    C ( 0.1, 1.0 ), 
    C ( 0.4, 0.6 ), 
    C ( 1.0, 3.0 ), 
    C ( 2.4, 2.0 ), 
    C (-0.4, 2.0 ), 
    C ( 1.1, 0.2 ), 
    C ( 0.8, 0.2 ), 
    C ( 0.3,-1.0 ) 
  };

  std::vector<C> Responses ( nDOFs * InputVars.size() / Dim );

  for ( auto i = 0; i < Responses.size(); i++ ) {

    Responses[2*i  ] = C ( 1.0, 0.0 ) / (
      InputVars[2*i  ] +
      InputVars[2*i+1]
    );

    Responses[2*i+1] = C ( 1.0, 0.0 ) / (
      InputVars[2*i  ]
    );

  }

  ASSERT_NO_THROW (
    Model.Train ( InputVars, Responses );
  );


  #ifdef DEBUG_RSMSD 

  auto numCoeffs = Model.NumCoeffs(); 
  auto denCoeffs = Model.DenCoeffs(); 

  R tolerance = 1e-4; 

  EXPECT_NEAR ( numCoeffs[0].real(), denCoeffs[1].real(), tolerance );
  EXPECT_NEAR ( numCoeffs[0].real(), denCoeffs[2].real(), tolerance );
  EXPECT_NEAR ( denCoeffs[1].real(), denCoeffs[2].real(), tolerance );

  EXPECT_NEAR ( numCoeffs[0].imag(), denCoeffs[1].imag(), tolerance );
  EXPECT_NEAR ( numCoeffs[0].imag(), denCoeffs[2].imag(), tolerance );
  EXPECT_NEAR ( denCoeffs[1].imag(), denCoeffs[2].imag(), tolerance );

  EXPECT_NEAR ( denCoeffs[0].real(), 0.0, tolerance );
  EXPECT_NEAR ( denCoeffs[0].imag(), 0.0, tolerance );

  EXPECT_NEAR ( numCoeffs[1].real(), denCoeffs[4].real(), tolerance );
  EXPECT_NEAR ( numCoeffs[1].imag(), denCoeffs[4].imag(), tolerance );

  EXPECT_NEAR ( denCoeffs[3].real(), 0.0, tolerance );
  EXPECT_NEAR ( denCoeffs[3].imag(), 0.0, tolerance );

  EXPECT_NEAR ( denCoeffs[5].real(), 0.0, tolerance );
  EXPECT_NEAR ( denCoeffs[5].imag(), 0.0, tolerance );

  #endif // DEBUG_RSMSD 


} // NIRPCE_Train, ComplexNDOF2_DIM2  


