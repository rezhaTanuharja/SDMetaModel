/**
 * @file substructuredModel.cpp
 *
 * @brief
 * Contains implementations of MCS using the UCB method.
 *
 * @author Rezha Adrian Tanuharja <rezhadr@outlook.com>
 * @date 2023-11-02
*/

#include "alias.hpp"
#include "materials.hpp"
#include "plates.hpp"
#include "solvers.hpp"
#include "substructuring.hpp"
#include "utilities.hpp"

#include <mpi.h>

int main ( int argc, char** argv ) {

  // ----- Initiate MPI environment -------------------------------------------

  MPI_Init ( &argc, &argv );

  int rank; 
  MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

  int size;
  MPI_Comm_size ( MPI_COMM_WORLD, &size );

  // ----- Define HUCB parameters ---------------------------------------------

  INT nSamples = 300;
  INT nOmega = 18;
  INT nModes[] = { 4, 16, 64 };

  FLOAT dOmega = 0.5;
  FLOAT omega0 = 9.0;

  // ----- Define plate's parameters ------------------------------------------

  // Substructured model consists of three components. The mesh parameters
  // of each component are such that the combined grid is equivalent to the
  // grid in the complete model.

  std::array<INT, 2> numEls;
  std::array<FLOAT, 2> lengths;
  std::array<FLOAT, 2> origins;
  std::vector<INT> boundaryDOFs;
  std::vector<INT> couplingDOFs;

  if ( rank == 0 ) {
    numEls = { 14, 20 };
    lengths = {  7.0 / 3.0, 4.0 };
    origins = {  8.0 / 3.0, 0.0 };
    boundaryDOFs = {      };
    couplingDOFs = { 0, 2 };
  }

  if ( rank == 1 ) {
    numEls = { 16, 20 };
    lengths = {  8.0 / 3.0, 4.0 };
    origins = {  0.0 / 3.0, 0.0 };
    boundaryDOFs = { 0    };
    couplingDOFs = {    2 };
  }

  if ( rank == 2 ) {
    numEls = { 18, 20 };
    lengths = {  9.0 / 3.0, 4.0 };
    origins = { 15.0 / 3.0, 0.0 };
    boundaryDOFs = {    2 };
    couplingDOFs = { 0    };
  }

  mFEM::RectangularPlate plate (
    numEls, lengths, origins, mFEM::gaussQuadrature
  );

  auto boundaryDofIds = plate.boundaryDofIds ( boundaryDOFs );
  auto couplingDofIds = plate.boundaryDofIds ( couplingDOFs );

  auto internalDofIds = plate.internalDofIds ( couplingDofIds );

  auto constraintDofIds = plate.boundaryDofIds ( { 0, 2 } );
  auto freeDofIds = plate.internalDofIds ( constraintDofIds );

  // ----- Define reference material ------------------------------------------

  mFEM::PlateMaterial steel;

  steel.density           = 7850.0;
  steel.shearCorrFactor   = 0.83;
  steel.shearModulus      = 79.3e9;
  steel.elasticityModulus = 200e9;
  steel.bendingStiffness  = 18315.0;
  steel.poissonRatio      = 0.3;
  steel.thickness         = 0.01;

  // ----- Define node id of unit load ----------------------------------------

  // The node is located inside component 2.

  INT nodeA = 3 * ( 2 * 2 + 1 ) * ( 2 * 20 + 1 ) - 3;

  // ----- Compute reference mass and stiffness matrices ----------------------

  mFEM::GlobalSystem referenceSystem = plate.assembleGlobalSystem ( 
    steel, [=](FLOAT, FLOAT) { return 0.0; }
  );

  mFEM::GlobalMatrix Mii_ref = mFEM::extract (
    std::get<0>( referenceSystem ),
    freeDofIds,
    freeDofIds
  );

  mFEM::GlobalMatrix Kii_ref = mFEM::extract (
    std::get<0>( referenceSystem ),
    freeDofIds,
    freeDofIds
  );

  // ----- Perform MCS using substructured model ------------------------------

  for ( INT i = 0; i < nOmega; i++ ) {

    FLOAT omega = omega0 + i * dOmega;

    for ( const INT mode : nModes ) {

      // ----- Compute reference normal modes ---------------------------------

      mFEM::MatrixXF Psi_free = mFEM::computeEigenvectors (
        Mii_ref, Kii_ref, omega, mode
      );

      mFEM::MatrixXF Psi_n = mFEM::zeroPadRows (
        Psi_free, freeDofIds, boundaryDofIds
      );

      for ( INT j = 0; j < nSamples; j++ ) {

        FLOAT mult1 = 1.0 + 0.2 * mFEM::normalVariable ( 0.0, 1.0, j + 1 );
        FLOAT mult2 = 1.0 + 0.2 * mFEM::normalVariable ( 0.0, 1.0, j + 7 );

        steel.density   = mult1 * 7850.0;
        steel.thickness = mult2 * 0.01;

        mFEM::GlobalSystem globalSystem = plate.assembleGlobalSystem ( 
          steel, [=](FLOAT, FLOAT) { return 0.0; }
        );

        mFEM::GlobalMatrix dynamicStiffness = mFEM::computeDynamicStiffness (
          std::get<0> ( globalSystem ),
          std::get<1> ( globalSystem ),
          omega
        );

        mFEM::GlobalVector unitLoad ( std::get<2>(globalSystem).size() );
        unitLoad = 0.0 * unitLoad;

        // Node A is inside component 2.

        if ( rank == 2 ) {
          unitLoad ( nodeA ) = -1.0 * omega * omega;
        }

        mFEM::applyHomogenousDirichletBC (
          dynamicStiffness,
          unitLoad,
          boundaryDofIds
        );

        // ----- Compute constraint modes -------------------------------------

        mFEM::GlobalMatrix Dii = mFEM::extract (
          dynamicStiffness, internalDofIds, internalDofIds 
        );

        mFEM::GlobalMatrix Dic = mFEM::extract (
          dynamicStiffness, internalDofIds, couplingDofIds 
        );

        Eigen::SparseLU<mFEM::GlobalMatrix> solver;

        solver.analyzePattern(Dii);
        solver.factorize(Dii);

        mFEM::MatrixXF Psi_c = solver.solve(-Dic);

        // ----- Perform CB transformation ------------------------------------

        auto subcomponent = mFEM::reduceComponent (
          dynamicStiffness, unitLoad,
          internalDofIds, couplingDofIds, Psi_n, Psi_c
        );

        // ----- Condense matrices and vectors on the boundaries --------------

        mFEM::MatrixXF S = std::get<0>( subcomponent );
        mFEM::VectorXF F = std::get<1>( subcomponent );

        std::vector<INT> globalSSizes ( size );
        INT localSSize = S.size();

        std::vector<INT> globalFSizes ( size );
        INT localFSize = F.size();

        INT nXb = 6 * numEls[1] + 3;

        mFEM::VectorXF U ( 2 * nXb );

        MPI_Allgather (
          &localSSize, 1, MPI_INT, globalSSizes.data(), 1, MPI_INT,
          MPI_COMM_WORLD
        );

        MPI_Allgather (
          &localFSize, 1, MPI_INT, globalFSizes.data(), 1, MPI_INT,
          MPI_COMM_WORLD
        );

        if ( rank != 0 ) {
          MPI_Send (
            S.data(), S.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD 
          );
          MPI_Send (
            F.data(), nXb, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD 
          );
        }

        if ( rank == 0 ) {
          
          std::vector<FLOAT> receivedS;
          std::vector<FLOAT> receivedF;

          receivedS.resize ( nXb * nXb );

          MPI_Recv (
            receivedS.data(), nXb * nXb, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD,
            MPI_STATUS_IGNORE
          );

          Eigen::Map<mFEM::MatrixXF> S1 ( 
            receivedS.data(), nXb, nXb
          );

          S.block ( 0, 0, nXb, nXb ) += S1;

          MPI_Recv (
            receivedS.data(), nXb * nXb, MPI_DOUBLE, 2, 0, MPI_COMM_WORLD,
            MPI_STATUS_IGNORE
          );

          Eigen::Map<mFEM::MatrixXF> S2 ( 
            receivedS.data(), nXb, nXb
          );

          S.block ( nXb, nXb, nXb, nXb ) += S2;

          receivedF.resize (nXb);

          MPI_Recv (
            receivedF.data(), nXb, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD,
            MPI_STATUS_IGNORE
          );

          Eigen::Map<mFEM::VectorXF> F1 ( 
            receivedF.data(), globalFSizes[1]
          );

          F.segment( 0, nXb ) += F1;

          MPI_Recv (
            receivedF.data(), nXb, MPI_DOUBLE, 2, 1, MPI_COMM_WORLD,
            MPI_STATUS_IGNORE
          );

          Eigen::Map<mFEM::VectorXF> F2 ( 
            receivedF.data(), globalFSizes[1]
          );

          F.segment( nXb, nXb ) += F2;

          U = S.inverse() * F;

        }

        MPI_Bcast (
          U.data(), 2*nXb, MPI_DOUBLE, 0, MPI_COMM_WORLD
        );

        mFEM::MatrixXF Dnn = std::get<2>( subcomponent );
        mFEM::MatrixXF Dnb = std::get<3>( subcomponent );
        mFEM::VectorXF Fn  = std::get<4>( subcomponent );

        Fn = 0.0 * Fn;

        mFEM::VectorXF xb;

        if ( rank == 0 ) {
          xb = U;
        }

        if ( rank == 1 ) {
          xb = U.segment(0,nXb);
        }

        if ( rank == 2 ) {
          xb = U.segment(nXb,nXb);
        }


        mFEM::VectorXF xi = Dnn.inverse() * ( Fn - Dnb * xb );

        mFEM::VectorXF xn = Psi_c * xb + Psi_n * xi;

        std::vector<FLOAT> X ( xn.size() + xb.size() );

        Eigen::Map<mFEM::VectorXF> x ( X.data(), xn.size() + xb.size() );

        x = mFEM::combine (
          xn, xb, internalDofIds, couplingDofIds
        );

        std::vector<FLOAT> disp ( 
          x.data(),
          x.data() + x.size() 
        );

        mFEM::saveVector (
          "HUCB/HUCB_" + std::to_string (rank) + "_" + std::to_string ( i ),
          disp 
        );

      } // loop through samples
    
    } // loop through number of modes

  } // loop through omega


  MPI_Finalize ();

} // main

