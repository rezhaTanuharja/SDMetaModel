/**
 * @file sparseRPCE.cpp
 *
 * @brief
 * Contains implementations of surrogate modelling techniques w/ sparse RPCE.
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

  FLOAT eps_0 = 1.0;

  INT nDOFs = 1;
  INT dimension = 2;

  INT nOmega = 18;
  INT nReps  = 40;

  std::vector<INT> seedOffsets = { 1001, 1007 };

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

    INT nPoint = trainData.size() / nDOFs;

    for ( INT j = 40; j < nPoint; j += 10 ) {
      for ( INT k = 0; k < nReps; k++ ) {

        std::vector<INT> pointID ( nPoint );
        std::iota ( pointID.begin(), pointID.end(), 0 );

        std::vector<INT> testID = surrogate::randomSubset ( pointID, j, k );

        std::vector<FLOAT> trainInputParameters;
        for ( INT p = 0; p < testID.size(); p++ ) {
          for ( INT q = 0; q < dimension; q++ ) {
            trainInputParameters.push_back (
              surrogate::normalVariable ( 
                0.0, 1.0, 
                testID[p] + seedOffsets[q] 
              )
            );
          }
        }

        std::vector<FLOAT> trainResponses;
        for ( INT p = 0; p < testID.size(); p++ ) {
          for ( INT q = 0; q < nDOFs; q++ ) {
            trainResponses.push_back (
              trainData[ p * nDOFs + q ]
            );
          }
        }

        model.trainSparse ( trainInputParameters, trainResponses, eps_0 );

        std::vector<INT> seeds;
        surrogate::loadVector (
          "MCS/seeds_" + std::to_string ( i ),
          seeds
        );

        std::vector<FLOAT> testInputParameters;
        for ( INT j = 0; j < seeds.size(); j++ ) {
          testInputParameters.push_back (
            surrogate::normalVariable ( 0.0, 1.0, seeds[j] + 1001 )
          );
          testInputParameters.push_back (
            surrogate::normalVariable ( 0.0, 1.0, seeds[j] + 1007 )
          );
        }

        std::vector<FLOAT> estimates = model.computeResponse ( 
          testInputParameters 
        );

        surrogate::saveVector (
          "sRPCE/sRPCE_" + node + "_64_" 
          + std::to_string ( i ) + "_"
          + std::to_string ( j ) + "_"
          + std::to_string ( k ),
          estimates
        );

      } // loop through repetitions (k)

    } // loop through number of samples (j)

  } // loop through omega (i)


} // main

