/**
 * @file materials.hpp
 *
 * @brief
 * Contains declarations of material properties struct.
 *
 * @author Rezha Adrian Tanuharja <rezhadr@outlook.com>
 * @date 2023-11-02
*/

#ifndef MATERIALS_DECLARATIONS
#define MATERIALS_DECLARATIONS

#ifndef ALIAS_DECLARATIONS
#include "alias.hpp"
#endif


namespace mFEM {


/**
 * @struct plateProperties
 *
 * @brief
 * Struct containing necessary material properties for plate elements.
*/
struct PlateMaterial {

  FLOAT density;
  FLOAT shearCorrFactor;
  FLOAT elasticityModulus;
  FLOAT shearModulus;
  FLOAT bendingStiffness;
  FLOAT poissonRatio;
  FLOAT thickness;

  PlateMaterial ();

}; // struct PlateMaterial


} // namespace mFEM


#endif // MATERIALS_DECLARATIONS

