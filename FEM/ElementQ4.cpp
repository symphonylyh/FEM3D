/**
 * @file ElementQ4.cpp
 * Implementation of ElementQ4 class.
 *
 * @author Haohang Huang
 * @date Feburary 10, 2018
 */

#include "ElementQ4.h"

// static member should be initialized outside the class body because it doesn't
// depend on any instance of that class and is initialized before any instance is
// created. Initialization here does not require the static member to be public.
ElementQ4::staticMembers ElementQ4::statics(4,4,4,2,2);

ElementQ4::ElementQ4()
{
}

ElementQ4::ElementQ4(const int & index, const std::vector<int> & nodeList, Node** const meshNode, Material* const material)
  : Element(index, nodeList, meshNode, material) // call the constructor of base class in the initializer list!
{
    if (!material->anisotropy) {
        modulusAtGaussPt = VectorXd::Constant(statics.shape->gaussianPt().size(), material_->modulus()); // 4 x 1 vector, the length depends on element type, so can only be initialized in derived class
    } else {
        modulusAtGaussPt = MatrixXd(statics.shape->gaussianPt().size(), 3); // 4 x 3 Matrix. Each column: horizontal modulus, vertical modulus, shear modulus
        modulusAtGaussPt.col(0) = VectorXd::Constant(statics.shape->gaussianPt().size(), material_->modulusR());
        modulusAtGaussPt.col(1) = VectorXd::Constant(statics.shape->gaussianPt().size(), material_->modulusZ());
        modulusAtGaussPt.col(2) = VectorXd::Constant(statics.shape->gaussianPt().size(), material_->modulusG());
    }
}

ElementQ4::~ElementQ4()
{
}

Shape* ElementQ4::shape() const
{
    return statics.shape;
}
