/**
 * @file ElementB20.cpp
 * Implementation of ElementB20 class.
 *
 * @author Haohang Huang
 * @date Feburary 13, 2018
 */

#include "ElementB20.h"

// static member should be initialized outside the class body because it doesn't
// depend on any instance of that class and is initialized before any instance is
// created. Initialization here does not require the static member to be public.
ElementB20::staticMembers ElementB20::statics(8,9,4,3,3);

ElementB20::ElementB20()
{
}

ElementB20::ElementB20(const int & index, const std::vector<int> & nodeList, Node** const meshNode, Material* const material)
    : Element(index, nodeList, meshNode, material) // call the constructor of base class in the initializer list!
{
    if (!material->anisotropy) {
        modulusAtGaussPt = VectorXd::Constant(statics.shape->gaussianPt().size(), material_->modulus()); // 9 x 1 vector, the length depends on element type, so can only be initialized in derived class
    } else {
        modulusAtGaussPt = MatrixXd(statics.shape->gaussianPt().size(), 3); // 9 x 3 Matrix. Each column: horizontal modulus, vertical modulus, shear modulus
        modulusAtGaussPt.col(0) = VectorXd::Constant(statics.shape->gaussianPt().size(), material_->modulusR());
        modulusAtGaussPt.col(1) = VectorXd::Constant(statics.shape->gaussianPt().size(), material_->modulusZ());
        modulusAtGaussPt.col(2) = VectorXd::Constant(statics.shape->gaussianPt().size(), material_->modulusG());
    }
}

ElementB20::~ElementB20()
{
}

Shape* ElementB20::shape() const
{
    return statics.shape;
}