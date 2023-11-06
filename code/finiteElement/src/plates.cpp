/**
 * @file plates.cpp
 *
 * @brief
 * Contains implementations of plate element class.
 *
 * @author Rezha Adrian Tanuharja <rezhadr@outlook.com>
 * @date 2023-11-02
*/

#ifndef PLATES_DECLARATIONS
#include "plates.hpp"
#endif


namespace mFEM {


RectangularPlate::RectangularPlate (
  std::array<INT,2> numberOfElements, 
  std::array<FLOAT,2> lengths, 
  std::array<FLOAT,2> origins, 
  IntegrationPointsProvider integrationPointsProvider 
) :
  numberOfElements_ ( numberOfElements ),
  lengths_ ( lengths ), 
  origins_ ( origins ), 
  integrationPointsProvider_ ( integrationPointsProvider )
{ 
  locationMaps_ = constructLocationMaps (); 
} 


const std::vector<INT> RectangularPlate::boundaryDofIds (
  const std::vector<INT> sides
) const {

  std::vector<INT> boundaryDofIds; 

  boundaryDofIds.reserve (
    12 * numberOfElements_[0] + 12 * numberOfElements_[1] + 12
  );

  for ( INT i = 0; i < sides.size(); i++ ) {

    if ( sides[i] == 0 ) {
      for ( INT j = 0; j < 2 * numberOfElements_[1] + 1; j++ ) {
        boundaryDofIds.push_back ( 3 * j     );
        boundaryDofIds.push_back ( 3 * j + 1 );
        boundaryDofIds.push_back ( 3 * j + 2 );
      }
      continue;
    }

    if ( sides[i] == 1 ) {
      for ( INT j = 0; j < 2 * numberOfElements_[0] + 1; j++ ) {
        boundaryDofIds.push_back (
          6 * j * numberOfElements_[1] + 3 * j
        );
        boundaryDofIds.push_back (
          6 * j * numberOfElements_[1] + 3 * j + 1
        );
        boundaryDofIds.push_back (
          6 * j * numberOfElements_[1] + 3 * j + 2
        );
      }
      continue;
    }

    if ( sides[i] == 2 ) {
      INT elementId = 
        2 * numberOfElements_[0] * ( 2 * numberOfElements_[1] + 1 );
      for ( INT j = 0; j < 2 * numberOfElements_[1] + 1; j++ ) {
        boundaryDofIds.push_back ( 3 * ( j + elementId )     );
        boundaryDofIds.push_back ( 3 * ( j + elementId ) + 1 );
        boundaryDofIds.push_back ( 3 * ( j + elementId ) + 2 );
      }
      continue;
    }

    if ( sides[i] == 3 ) {
      INT elementId = 2 * numberOfElements_[1];
      for ( INT j = 0; j < 2 * numberOfElements_[0] + 1; j++ ) {
        boundaryDofIds.push_back (
          3 * ( elementId + 2 * j * numberOfElements_[1] + j )
        );
        boundaryDofIds.push_back (
          3 * ( elementId + 2 * j * numberOfElements_[1] + j ) + 1
        );
        boundaryDofIds.push_back (
          3 * ( elementId + 2 * j * numberOfElements_[1] + j ) + 2
        );
      }
      continue;
    }

  }
  
  std::sort ( boundaryDofIds.begin(), boundaryDofIds.end() );
  boundaryDofIds.erase (
    std::unique ( boundaryDofIds.begin(), boundaryDofIds.end() ),
    boundaryDofIds.end()
  );

  return boundaryDofIds;

} // RectangularPlate::boundaryDofIds 


const std::vector<INT> RectangularPlate::internalDofIds (
  const std::vector<INT>& boundaryDofIds 
) const {

  INT nDOFs = 3 * 
    ( 2 * numberOfElements_[0] + 1 ) * 
    ( 2 * numberOfElements_[1] + 1 );

  std::vector<INT> internalDofIds ( nDOFs );
  std::iota ( internalDofIds.begin(), internalDofIds.end(), 0 );

  internalDofIds.erase (
    std::remove_if (
      internalDofIds.begin(), 
      internalDofIds.end(), 
      [&] ( int num ) {
        return std::binary_search (
          boundaryDofIds.begin(), 
          boundaryDofIds.end(), 
          num
        );
      }
    ), 
    internalDofIds.end()
  );

  return internalDofIds;

} // RectangularPlate::internalDofIds


const std::array<const std::vector<FLOAT>,2> RectangularPlate::nodeCoordinates () 
const {

  INT nX = 2 * numberOfElements_[0] + 1;
  INT nY = 2 * numberOfElements_[1] + 1;
  INT nDOFs = nX * nY;

  std::vector<FLOAT> X;
  X.reserve ( nDOFs );

  std::vector<FLOAT> Y;
  Y.reserve ( nDOFs );

  for ( INT i = 0; i < nX; i++ ) {
    for ( INT j = 0; j < nY; j++ ) {
      X.push_back (
        origins_[0] + 0.5 * i * lengths_[0] / numberOfElements_[0]
      );
    }
  }

  for ( INT i = 0; i < nX; i++ ) {
    for ( INT j = 0; j < nY; j++ ) {
      Y.push_back (
        origins_[1] + 0.5 * j * lengths_[1] / numberOfElements_[1]
      );
    }
  }

  return { X, Y };

} // RectangularPlate::nodeCoordinates 


const LocationMaps RectangularPlate::locationMaps () const {
  return locationMaps_;
}


GlobalSystem RectangularPlate::assembleGlobalSystem (
  const PlateMaterial& material, 
  const SpatialFunction& force
) const {

  INT nDOFs = 3 * 
    ( 2 * numberOfElements_[0] + 1 ) *
    ( 2 * numberOfElements_[1] + 1 );

  GlobalMatrix globalMass; 
  GlobalMatrix globalStiffness;
  GlobalVector globalLoad = GlobalVector::Zero ( nDOFs );

  globalMass.resize ( nDOFs, nDOFs );
  globalStiffness.resize ( nDOFs, nDOFs );

  std::vector<Eigen::Triplet<FLOAT>> mass;
  std::vector<Eigen::Triplet<FLOAT>> stiffness;

  INT nY = numberOfElements_[1];

  for ( INT iEl = 0; iEl < numberOfElements_[0]; iEl++ ) {
    for ( INT jEl = 0; jEl < numberOfElements_[1]; jEl++ ) {

      INT elIndex = iEl * nY + jEl;

      ElementSystem elementSystem = integrateElementSystem (
        { iEl, jEl }, material, force
      );

      for ( INT i = 0; i < 27; i++ ) {
        for ( INT j = 0; j < 27; j++ ) {

          FLOAT massValue = std::get<0>(elementSystem)[i*27+j];

          mass.emplace_back (
            locationMaps_[elIndex][i],
            locationMaps_[elIndex][j],
            massValue
          );

          FLOAT stiffnessValue = std::get<1>(elementSystem)[i*27+j];

          stiffness.emplace_back (
            locationMaps_[elIndex][i],
            locationMaps_[elIndex][j],
            stiffnessValue
          );

        }
      }

      for ( INT i = 0; i < 27; i++ ) {
        globalLoad ( locationMaps_[elIndex][i] ) += 
          std::get<2>(elementSystem)[i];
      }

    }
  }

  globalMass.setFromTriplets ( mass.begin(), mass.end() );
  globalStiffness.setFromTriplets ( stiffness.begin(), stiffness.end() );

  return { globalMass, globalStiffness, globalLoad };

} // RectangularPlate::assembleGlobalSystem 


LocationMaps RectangularPlate::constructLocationMaps () {

  LocationMaps locationMaps;

  INT nDOFs1 = 6 * numberOfElements_[1] + 3; 

  for ( INT iElement = 0; iElement < numberOfElements_[0]; iElement++ ) {
    for ( INT jElement = 0; jElement < numberOfElements_[1]; jElement++ ) {

      INT firstId = iElement * 2 * nDOFs1 + 6 * jElement;

      LocationMap locationMap;

      for ( INT i = 0; i < 3; i++ ) {
        for ( INT j = 0; j < 9; j++ ) {
          locationMap.push_back ( firstId + j + i * nDOFs1 );
        }
      } // loop over nodes in biquadratic element

      locationMaps.push_back ( locationMap );

    }
  } // loop over all elements in 2D mesh

  return locationMaps;

} // RectangularPlate::constructLocationMaps 


const std::vector<FLOAT> RectangularPlate::evaluateActiveBasis (
  const std::array<FLOAT,2> localCoordinates, 
  const std::array<INT,2> diffOrders
) const {

  std::vector<FLOAT> basisFunctions; 

  auto N = [=]( const INT dimension, const INT functionIndex ) {
    return quadraticBasis ( 
      localCoordinates[dimension],
      functionIndex,
      diffOrders[dimension]
    );
  };

  for ( INT i = 0; i < 3; i++ ) {
    for ( INT j = 0; j < 3; j++ ) {
      basisFunctions.push_back (
        N (0,i) * N (1,j)
      );
    }
  }

  return basisFunctions;

} // RectangularPlate::evaluateActiveBasis 


const std::array<FLOAT,2> RectangularPlate::mapToGlobalCoordinates (
  const std::array<FLOAT,2> localCoordinates, 
  const std::array<INT,2> elementIndices
) const {
  
  auto coordinates = [=] ( const INT i ) {
    return origins_[i] + lengths_[i] / numberOfElements_[i] * 
      ( 0.5 * ( 1 + localCoordinates[i] ) + elementIndices[i] );
  };

  return { coordinates(0), coordinates(1) };

} // mFEM::RectangularPlate::mapToGlobalCoordinates 


ElementSystem RectangularPlate::integrateElementSystem (
  const std::array<INT,2> elementIndices,
  const PlateMaterial& material, 
  const SpatialFunction& force
) const {

  INT nDOFs = 27; 

  std::vector<FLOAT> massMatrix      ( nDOFs * nDOFs, 0.0 );
  std::vector<FLOAT> stiffnessMatrix ( nDOFs * nDOFs, 0.0 );
  std::vector<FLOAT> forceVector     (         nDOFs, 0.0 );

  FLOAT alpha = material.shearCorrFactor;
  FLOAT G     = material.shearModulus;
  FLOAT E     = material.elasticityModulus;
  FLOAT t     = material.thickness;
  FLOAT nu    = material.poissonRatio; 
  FLOAT rho   = material.density;
  FLOAT D     = material.bendingStiffness; 
  // FLOAT D = E * t * t * t * t / (12.0 * (1.0 - nu * nu) );

  Eigen::Matrix<FLOAT,3,3> I = Eigen::Matrix<FLOAT,3,3>::Zero();

  I (0,0) = rho * t;
  I (1,1) = rho * t * t * t / 12.0;
  I (2,2) = rho * t * t * t / 12.0;

  Eigen::Matrix<FLOAT,5,5> C = Eigen::Matrix<FLOAT,5,5>::Zero();

  C (0,0) = alpha * G * t;
  C (1,1) = alpha * G * t;
  C (2,2) = D;
  C (2,3) = nu * D;
  C (3,2) = nu * D;
  C (3,3) = D;
  C (4,4) = 0.5 * D * ( 1.0 - nu );

  FLOAT Jx = 0.5 * lengths_[0] / numberOfElements_[0];
  FLOAT Jy = 0.5 * lengths_[1] / numberOfElements_[1];
  FLOAT detJ = Jx * Jy;

  Eigen::Map<ElementMatrix> M (      massMatrix.data(), 27, 27 );
  Eigen::Map<ElementMatrix> K ( stiffnessMatrix.data(), 27, 27 );
  Eigen::Map<ElementVector> F (     forceVector.data(),     27 );

  IntegrationPoints integrationPoints = integrationPointsProvider_ (3);

  for ( INT i = 0; i < integrationPoints[0].size(); i++ ) {
    for ( INT j = 0; j < integrationPoints[0].size(); j++ ) {

      FLOAT weight = integrationPoints[1][i] * integrationPoints[1][j];

      FLOAT xi = integrationPoints[0][i];
      FLOAT et = integrationPoints[0][j];

      auto  shapeFunctions    = evaluateActiveBasis ( {xi,et}, {0,0} );
      auto DshapeFunctionsDxi = evaluateActiveBasis ( {xi,et}, {1,0} );
      auto DshapeFunctionsDet = evaluateActiveBasis ( {xi,et}, {0,1} );

      Eigen::Matrix<FLOAT,3,27> N = Eigen::Matrix<FLOAT,3,27>::Zero();

      for ( INT k = 0; k < shapeFunctions.size(); k++ ) {
        N ( 0, 3*k     ) = shapeFunctions[k];
        N ( 1, 3*k + 1 ) = shapeFunctions[k];
        N ( 2, 3*k + 2 ) = shapeFunctions[k];
      }

      Eigen::Matrix<FLOAT,5,27> B = Eigen::Matrix<FLOAT,5,27>::Zero();

      for ( INT k = 0; k < shapeFunctions.size(); k++ ) {
        B ( 0, 3*k     ) =  DshapeFunctionsDxi[k] / Jx;
        B ( 1, 3*k     ) =  DshapeFunctionsDet[k] / Jy;
        B ( 1, 3*k + 1 ) = - shapeFunctions   [k]     ;
        B ( 3, 3*k + 1 ) = -DshapeFunctionsDet[k] / Jy;
        B ( 4, 3*k + 1 ) = -DshapeFunctionsDxi[k] / Jx;
        B ( 0, 3*k + 2 ) =   shapeFunctions   [k]     ;
        B ( 2, 3*k + 2 ) =  DshapeFunctionsDxi[k] / Jx;
        B ( 4, 3*k + 2 ) =  DshapeFunctionsDet[k] / Jy;
      }

      auto globalCoordinates = mapToGlobalCoordinates (
        { xi, et }, elementIndices
      );

      Eigen::Vector<FLOAT,3> Load = Eigen::Vector<FLOAT,3>::Zero();
      Load (0) = force ( globalCoordinates[0], globalCoordinates[1] );
      Load (1) = 0.0;
      Load (2) = 0.0;

      M += detJ * weight * N.transpose() * I * N;
      K += detJ * weight * B.transpose() * C * B;
      F += detJ * weight * N.transpose() * Load;

    }
  }

  return { massMatrix, stiffnessMatrix, forceVector };

} // RectangularPlate::integrateElementSystem 


} // namespace mFEM

