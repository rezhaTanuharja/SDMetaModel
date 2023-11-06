/**
 * @file plates.hpp
 *
 * @brief
 * Contains declarations of plate element class.
 *
 * @author Rezha Adrian Tanuharja <rezhadr@outlook.com>
 * @date 2023-11-02
*/

#ifndef PLATES_DECLARATIONS
#define PLATES_DECLARATIONS

#ifndef ALIAS_DECLARATIONS
#include "alias.hpp"
#endif

#ifndef BASIS_FUNCTIONS_DECLARATIONS
#include "basisFunctions.hpp"
#endif

#ifndef MATERIALS_DECLARATIONS
#include "materials.hpp"
#endif

#include <numeric>


namespace mFEM {


/**
 * @class rectangularPlate
 *
 * @brief
 * Shear deformable plate with regular biquadratic mesh.
*/
class RectangularPlate {

  std::array<INT, 2> numberOfElements_;
  std::array<FLOAT, 2> lengths_;
  std::array<FLOAT, 2> origins_;

  LocationMaps locationMaps_;

  IntegrationPointsProvider integrationPointsProvider_;

public: // class RectangularPlate


  /**
    * @brief 
    * Construct shear deformable plate with regular biquadratic mesh.
    *
    * @param numberOfElements           number of elements in x and y direction
    * @param lengths                    length of plate in x and y direction
    * @param origins                    coordinates of bottom left corner
    * @param integrationPointsProvider  function to generate integration points
    */
  RectangularPlate (
    std::array<INT,2> numberOfElements, 
    std::array<FLOAT,2> lengths, 
    std::array<FLOAT,2> origins, 
    IntegrationPointsProvider integrationPointsProvider 
  );


  /**
    * @brief 
    * Generate Dof ids along plate's side(s). 
    * @param sides vector of sides' numbers: 0(W), 1(S), 2(E), and/or 3(N)
    * @return vector of unique Dof ids in ascending order
    */
  const std::vector<INT> boundaryDofIds (
    const std::vector<INT> sides
  ) const;


  /**
    * @brief 
    * Generate Dof ids excluding the boundary Dof ids.
    * @param boundaryDofIds vector of boundary Dof ids
    * @return vector of unique Dof ids in ascending order 
    */
  const std::vector<INT> internalDofIds (
    const std::vector<INT>& boundaryDofIds
  ) const;


  /**
    * @brief 
    * Generate node coordinates. 
    * @return vector of x- and y- coordinates in column major order
    */
  const std::array<const std::vector<FLOAT>,2> nodeCoordinates () const;


  /**
    * @brief 
    * Retrieve mappings between element Dof ids and global Dof ids.
    * @return vector of element mappings
    */
  const LocationMaps locationMaps () const;


  /**
    * @brief 
    * Compute global mass matrix, stiffness matrix, and force vector.
    * @param properties   plate properties
    * @param force        force as function of spatial location
    * @return sparse mass and stiffness matrix and dense load vector
    */
  GlobalSystem assembleGlobalSystem (
    const PlateMaterial& material, 
    const SpatialFunction& force
  ) const;


private: // class RectangularPlate


  /**
    * @private 
    * Construct mapping between local and global DoFs.
    */
  LocationMaps constructLocationMaps ();


  /**
    * @private 
    * Evaluate basis functions or their derivatives at given local coordinates.
    */
  const std::vector<FLOAT> evaluateActiveBasis (
    const std::array<FLOAT,2> localCoordinates, 
    const std::array<INT,2> diffOrders
  ) const;


  /**
    * @private 
    * Map local element coordinates to global plate coordinates.
    */
  const std::array<FLOAT,2> mapToGlobalCoordinates (
    const std::array<FLOAT,2> localCoordinates, 
    const std::array<INT,2> elementIndices
  ) const;


  /**
    * @private 
    * Compute element mass matrix, stiffness matrix, and load vector.
    */
  ElementSystem integrateElementSystem (
    const std::array<INT,2> elementIndices,
    const PlateMaterial& material, 
    const SpatialFunction& force
  ) const;


}; // class RectangularPlate


} // namespace mFEM


#endif // PLATES_DECLARATIONS

