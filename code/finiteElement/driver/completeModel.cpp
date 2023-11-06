/**
 * @file completeModel.cpp
 *
 * @brief
 * Contains implementations of direct MCS using the complete plate model.
 *
 * @author Rezha Adrian Tanuharja <rezhadr@outlook.com>
 * @date 2023-11-02
*/

#include "alias.hpp"
#include "materials.hpp"
#include "plates.hpp"
#include "solvers.hpp"
#include "utilities.hpp"

int main ( int argc, char** argv ) {

  // ----- Define direct MCS parameters ---------------------------------------

  INT nSamples = 300;
  INT nOmega = 18;

  FLOAT dOmega = 0.5;
  FLOAT omega0 = 9.0;

  // ----- Define iterative solver's parameters -------------------------------

  INT maxIter = 250000;
  FLOAT tolerance = 1e-6;

  // ----- Define plate's parameters ------------------------------------------

  std::array<INT, 2> numEls = { 48, 20 };
  std::array<FLOAT, 2> lengths = { 8.0, 4.0 };
  std::array<FLOAT, 2> origins = { 0.0, 0.0 };

  mFEM::RectangularPlate plate (
    numEls, lengths, origins, mFEM::gaussQuadrature
  );

  std::vector<INT> boundaryDOFs = { 0, 2 };

  auto boundaryDofIds = plate.boundaryDofIds ( boundaryDOFs );

  // ----- Define material ----------------------------------------------------

  mFEM::PlateMaterial steel;

  steel.density           = 7850.0;
  steel.shearCorrFactor   = 0.83;
  steel.shearModulus      = 79.3e9;
  steel.elasticityModulus = 200e9;
  steel.bendingStiffness  = 18315.0;
  steel.poissonRatio      = 0.3;
  steel.thickness         = 0.01;

  // ----- Define node id of unit load ----------------------------------------

  INT nodeA = 3 * ( 2 * 32 + 1 ) * ( 2 * 20 + 1 ) - 3;

  // ----- Perform direct MCS -------------------------------------------------

  for ( INT i = 0; i < nOmega; i++ ) {

    FLOAT omega = omega0 + i * dOmega;

    std::vector<INT> converged;

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
      unitLoad ( nodeA ) = -1.0 * omega * omega;

      mFEM::applyHomogenousDirichletBC (
        dynamicStiffness,
        unitLoad,
        boundaryDofIds
      );

      auto result = mFEM::solveDynamicSystem (
        dynamicStiffness,
        unitLoad,
        maxIter,
        tolerance
      );

      if ( result.second ) {

        std::vector<FLOAT> DOFs (
          result.first.data(),
          result.first.data() + result.first.size()
        );


        // Append solution to MCS_(omega)

        mFEM::saveVector (
          "MCS/MCS_" + std::to_string (i),
          DOFs 
        );

        converged.push_back ( j );

      }

    } // loop through samples


    // Save sample indices that converges

    mFEM::saveVector (
      "MCS/seeds_" + std::to_string ( i ),
      converged 
    );

  } // loop through omega

} // main
