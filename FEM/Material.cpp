/**
 * @file Material.cpp
 * Implementation of the Material class.
 *
 * @author Haohang Huang
 * @date Apr 27, 2018
 */

#include "Material.h"

Material::Material(const bool & Anisotropy, const bool & Nonlinearity, const bool & NoTension)
  : anisotropy(Anisotropy), nonlinearity(Nonlinearity), noTension(NoTension), E_(MatrixXd::Zero(4,4)), thermalStrain_(VectorXd::Zero(4))
{
}

Material::~Material()
{
}

const double & Material::modulus() const
{
    return M_;
}

const double & Material::modulusR() const
{
    return Mr_;
}

const double & Material::modulusZ() const
{
    return Mz_;
}

const double & Material::modulusG() const
{
    return G_;
}

void Material::adjustModulus(const double & ratio)
{
    (void)ratio; // silence warning
    return;
}

const MatrixXd & Material::EMatrix() const
{
    return E_;
}

VectorXd Material::stressDependentModulus(const VectorXd & stress) const
{
    (void)stress; // silence warning
    return VectorXd::Zero(3); // to silent warning
}

MatrixXd Material::EMatrix(const VectorXd & modulus) const
{
    (void)modulus; // silence warning
    return MatrixXd::Zero(4,4); // to silent warning
}

const Vector2d & Material::bodyForce() const
{
    return bodyForce_;
}

void Material::setBodyForce(const Vector2d & force)
{
    bodyForce_ = force;
}

const VectorXd & Material::thermalStrain() const
{
    return thermalStrain_;
}

void Material::setThermalStrain(const VectorXd & thermalStrain)
{
    thermalStrain_ = thermalStrain;
}
