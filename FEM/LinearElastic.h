/**
 * @file LinearElastic.h
 * Derived class from Material for linear elastic materials.
 *
 * @author Haohang Huang
 * @date May 19, 2018
 */

#ifndef LinearElastic_h
#define LinearElastic_h

#include "Material.h"

class LinearElastic : public Material
{
  public:
    /** See the documentation of base class Material. */
    LinearElastic(const bool & anisotropy, const bool & nonlinearity, const bool & noTension, const std::vector<double> & properties);
    ~LinearElastic();

    void adjustModulus(const double & ratio); // Used in back analysis scheme

};

#endif /* LinearElastic_h */
