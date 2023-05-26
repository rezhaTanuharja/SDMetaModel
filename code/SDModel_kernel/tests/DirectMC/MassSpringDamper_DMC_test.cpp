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

#include "MassSpringDamper_DMC.hpp" 
#include <gtest/gtest.h> 


TEST ( MSD_DMC, StaticUndamped_1DOF ) {

  typedef size_t Z; 
  typedef Stochastic::R R; 
  typedef Stochastic::C C; 
  
  Z nDOFs = 1; 
  Z Dim   = 2; 

  R Omega = 0.0; 

  Stochastic::DirectMC::MassSpringDamper Model ( nDOFs, Dim );

  std::vector<R> MassCoeffs {
    1.0, 0.0 
  };

  std::vector<R> SpringCoeffs {
    1.0, 0.5 
  };

  std::vector<C> ForceCoeffs {
    C ( 1.0, 0.0 ), 
    C ( 0.0, 0.0 )
  };

  Model.SetMassCoeffs ( MassCoeffs ); 
  Model.SetSpringCoeffs ( SpringCoeffs ); 
  Model.SetForceCoeffs ( ForceCoeffs ); 

  std::vector<R> RandomBasis {
    1.0, 2.0 
  }; 

  std::vector<C> expected {
    C ( 0.5, 0.0 )
  };

  auto result = Model.ComputeResponses ( RandomBasis, Omega );

  ASSERT_EQ ( result.size(), expected.size() ); 

  for ( auto i = 0; i < result.size(); i++ ) {
    EXPECT_EQ ( result[i].real(), expected[i].real() );
    EXPECT_EQ ( result[i].imag(), expected[i].imag() );
  }

} // MSD_DMC, StaticUndamped_1DOF 


TEST ( MSD_DMC, StaticDamped_1DOF ) {

  typedef size_t Z; 
  typedef Stochastic::R R; 
  typedef Stochastic::C C; 
  
  Z nDOFs = 1; 
  Z Dim   = 2; 

  R DampingRatio = 0.75; 
  R Omega = 0.0; 

  Stochastic::DirectMC::MassSpringDamper Model ( nDOFs, Dim );

  std::vector<R> MassCoeffs {
    1.0, 0.0 
  };

  std::vector<R> SpringCoeffs {
    1.0, 0.5 
  };

  std::vector<C> ForceCoeffs {
    C ( 1.0, 0.0 ), 
    C ( 0.0, 0.0 )
  };

  Model.SetMassCoeffs ( MassCoeffs ); 
  Model.SetSpringCoeffs ( SpringCoeffs ); 
  Model.SetForceCoeffs ( ForceCoeffs ); 

  Model.SetDampingRatio ( DampingRatio ); 

  std::vector<R> RandomBasis {
    3.0, 2.0 
  }; 

  std::vector<C> expected {
    C ( 0.75, 0.0 )
  };

  std::vector<C> result = Model.ComputeResponses ( RandomBasis, Omega );

  ASSERT_EQ ( result.size(), expected.size() ); 

  for ( auto i = 0; i < result.size(); i++ ) {
    EXPECT_FLOAT_EQ ( result[i].real(), expected[i].real() );
    EXPECT_FLOAT_EQ ( result[i].imag(), expected[i].imag() );
  }

} // MSD_DMC, StaticDamped_1DOF  


TEST ( MSD_DMC, DynamicUndamped_1DOF ) {

  typedef size_t Z; 
  typedef Stochastic::R R; 
  typedef Stochastic::C C; 
  
  Z nDOFs = 1; 
  Z Dim   = 2; 

  R DampingRatio = 0.0; 
  R Omega = 1.0; 

  Stochastic::DirectMC::MassSpringDamper Model ( nDOFs, Dim );

  std::vector<R> MassCoeffs {
    1.0, 0.0 
  };

  std::vector<R> SpringCoeffs {
    1.0, 0.5 
  };

  std::vector<C> ForceCoeffs {
    C ( 1.0, 0.0 ), 
    C ( 0.0, 0.0 )
  };

  Model.SetMassCoeffs ( MassCoeffs ); 
  Model.SetSpringCoeffs ( SpringCoeffs ); 
  Model.SetForceCoeffs ( ForceCoeffs ); 

  Model.SetDampingRatio ( DampingRatio ); 

  std::vector<R> RandomBasis {
    2.0, 2.0 
  }; 

  std::vector<C> expected {
    C ( 2.0, 0.0 )
  };

  std::vector<C> result = Model.ComputeResponses ( RandomBasis, Omega );

  ASSERT_EQ ( result.size(), expected.size() ); 

  for ( auto i = 0; i < result.size(); i++ ) {
    EXPECT_FLOAT_EQ ( result[i].real(), expected[i].real() );
    EXPECT_FLOAT_EQ ( result[i].imag(), expected[i].imag() );
  }

} // MSD_DMC, DynamicUndamped_1DOF  


TEST ( MSD_DMC, DynamicDamped_1DOF ) {

  typedef size_t Z; 
  typedef Stochastic::R R; 
  typedef Stochastic::C C; 
  
  Z nDOFs = 1; 
  Z Dim   = 2; 

  R DampingRatio = 0.3; 
  R Omega = 1.0; 

  Stochastic::DirectMC::MassSpringDamper Model ( nDOFs, Dim );

  std::vector<R> MassCoeffs {
    1.0, 0.0 
  };

  std::vector<R> SpringCoeffs {
    1.0, 0.5 
  };

  std::vector<C> ForceCoeffs {
    C ( 1.0, 0.0 ), 
    C ( 0.0, 0.0 )
  };

  Model.SetMassCoeffs ( MassCoeffs ); 
  Model.SetSpringCoeffs ( SpringCoeffs ); 
  Model.SetForceCoeffs ( ForceCoeffs ); 

  Model.SetDampingRatio ( DampingRatio ); 

  std::vector<R> RandomBasis {
    1.0, 8.0 
  }; 

  std::vector<C> expected {
    C ( 0.16, -0.12 )
  };

  std::vector<C> result = Model.ComputeResponses ( RandomBasis, Omega );

  ASSERT_EQ ( result.size(), expected.size() ); 

  for ( auto i = 0; i < result.size(); i++ ) {
    EXPECT_FLOAT_EQ ( result[i].real(), expected[i].real() );
    EXPECT_FLOAT_EQ ( result[i].imag(), expected[i].imag() );
  }

} // MSD_DMC, DynamicDamped_1DOF  


#ifdef DEBUG_RSMSD 


TEST ( MSD_DMC, ProperMassCoeffs ) {
  
  typedef size_t Z; 
  typedef Stochastic::R R; 
  
  Z nDOFs = 2; 
  Z Dim   = 2; 
  
  Stochastic::DirectMC::MassSpringDamper Model ( nDOFs, Dim );
  
  std::vector<R> MassCoeffs {
    1.0, 2.0, 3.0, 4.0 
  };

  Model.SetMassCoeffs ( MassCoeffs );

  auto result = Model.MassCoeffs(); 

  ASSERT_EQ ( MassCoeffs.size(), result.size() );

  for ( auto i = 0; i < result.size(); i++ ) {
    EXPECT_EQ ( MassCoeffs[i], result[i] );
  }
  
} // MSD_DMC, ProperMassCoeffs 


TEST ( MSD_DMC, ImproperMassCoeffs ) {
  
  typedef size_t Z; 
  typedef Stochastic::R R; 
  
  Z nDOFs = 2; 
  Z Dim   = 2; 
  
  Stochastic::DirectMC::MassSpringDamper Model ( nDOFs, Dim );
  
  std::vector<R> MassCoeffs {
    1.0, 2.0, 3.0 
  };

  bool exception_thrown = false; 

  try {

    Model.SetMassCoeffs ( MassCoeffs );

  } catch ( const std::exception& e ) {

    EXPECT_STREQ (
      "MSD_DMC::SetMassCoeffs: incorrect size", 
      e.what()
    );

    exception_thrown = true; 

  }

  EXPECT_TRUE ( exception_thrown );
  
} // MSD_DMC, ImproperMassCoeffs 


TEST ( MSD_DMC, ProperSpringCoeffs ) {
  
  typedef size_t Z; 
  typedef Stochastic::R R; 
  
  Z nDOFs = 2; 
  Z Dim   = 2; 
  
  Stochastic::DirectMC::MassSpringDamper Model ( nDOFs, Dim );
  
  std::vector<R> SpringCoeffs {
    1.0, 2.0, 3.0, 4.0 
  };

  Model.SetSpringCoeffs ( SpringCoeffs );

  auto result = Model.SpringCoeffs(); 

  ASSERT_EQ ( SpringCoeffs.size(), result.size() );

  for ( auto i = 0; i < result.size(); i++ ) {
    EXPECT_EQ ( SpringCoeffs[i], result[i] );
  }
  
} // MSD_DMC, ProperSpringCoeffs 


TEST ( MSD_DMC, ImproperSpringCoeffs ) {
  
  typedef size_t Z; 
  typedef Stochastic::R R; 
  
  Z nDOFs = 2; 
  Z Dim   = 2; 
  
  Stochastic::DirectMC::MassSpringDamper Model ( nDOFs, Dim );
  
  std::vector<R> SpringCoeffs {
    1.0, 2.0, 3.0 
  };

  bool exception_thrown = false; 

  try {

    Model.SetSpringCoeffs ( SpringCoeffs );

  } catch ( const std::exception& e ) {

    EXPECT_STREQ (
      "MSD_DMC::SetSpringCoeffs: incorrect size", 
      e.what()
    );

    exception_thrown = true; 

  }

  EXPECT_TRUE ( exception_thrown );
  
  
} // MSD_DMC, ImproperSpringCoeffs 


TEST ( MSD_DMC, NegDampingRatio ) {
  
  typedef size_t Z; 
  typedef Stochastic::R R; 
  
  Z nDOFs = 2; 
  Z Dim   = 2; 
  
  Stochastic::DirectMC::MassSpringDamper Model ( nDOFs, Dim );

  R DampingRatio = -0.5;
  
  bool exception_thrown = false; 

  try {

    Model.SetDampingRatio ( DampingRatio );

  } catch ( const std::exception& e ) {

    EXPECT_STREQ (
      "MSD_DMC::SetDampingRatio: negative DampingRatio", 
      e.what()
    );

    exception_thrown = true; 

  }

  EXPECT_TRUE ( exception_thrown );
  
  
} // MSD_DMC, NegDampingRatio  


#endif 

