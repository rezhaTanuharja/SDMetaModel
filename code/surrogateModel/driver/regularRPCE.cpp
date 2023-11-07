/**
 * @file regularRPCE.cpp
 *
 * @brief
 * Contains implementations of surrogate modelling techniques w/ regular RPCE.
 *
 * @author Rezha Adrian Tanuharja <rezhadr@outlook.com>
 * @date 2023-11-02
*/

#include "basisFunctions.hpp"
#include "models.hpp"
#include "utilities.hpp"


int main () {

  // ----- Set number of threads for OpenMP -----------------------------------

  omp_set_num_threads ( 8 );

  // ----- Set RPCE parameters ------------------------------------------------

  INT nDOFs = 1;
  INT dimension = 2;

  INT nOmega = 18;

  std::string node = "A";

  std::vector<INT> numIndices = surrogate::generateIndices ( dimension, 5 );
  std::vector<INT> denIndices = surrogate::generateIndices ( dimension, 6 );

  surrogate::RPCE model ( nDOFs, dimension );

  model.setNumIndices ( numIndices );
  model.setDenIndices ( denIndices );

  for ( INT i = 0; i < nOmega; i++ ) {

    std::vector<FLOAT> trainData;
    surrogate::loadVector (
      "HUCB/HUCB_" + node + "_64_" + std::to_string ( i ),
      trainData
    );

    std::vector<INT> seeds;
    surrogate::loadVector (
      "MCS/seeds_" + std::to_string ( i ),
      seeds
    );

    std::vector<FLOAT> inputParameters;
    for ( INT j = 0; j < seeds.size(); j++ ) {
      inputParameters.push_back (
        surrogate::normalVariable ( 0.0, 1.0, seeds[j] + 1 )
      );
      inputParameters.push_back (
        surrogate::normalVariable ( 0.0, 1.0, seeds[j] + 7 )
      );
    }

    model.train ( inputParameters, trainData );

    std::vector<FLOAT> estimates = model.computeResponse ( inputParameters );
    surrogate::saveVector (
      "rRPCE/rRPCE_" + node + "_64_" + std::to_string ( i ),
      estimates
    );

  } // loop through omega

} // main

