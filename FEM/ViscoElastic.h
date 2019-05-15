/**
 * @file ViscoElastic.h
 * Derived class from Material for viscoelastic materials.
 *
 * @author Jiayi Luo
 * @date May 19, 2018
 */

#ifndef ViscoElastic_h
#define ViscoElastic_h

#include "Material.h"

class ViscoElastic : public Material
{
  public:
    /** See the documentation of base class Material. */
    ViscoElastic(const bool & anisotropy, const bool & nonlinearity, const bool & noTension, const std::vector<double> & properties);
    ~ViscoElastic();
};

#endif /* ViscoElastic_h */
