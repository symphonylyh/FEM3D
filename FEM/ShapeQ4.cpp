/**
 * @file ShapeQ4.cpp
 * Implementation of ShapeQ8 class.
 *
 * @author Haohang Huang
 * @date Feburary 10, 2018
 */

#include "ShapeQ4.h"

ShapeQ4::ShapeQ4(const int & nodes, const int & gaussians, const int & edges, const int & edgeNodes, const int & edgeGaussians) :
    Shape(nodes, gaussians, edges, edgeNodes, edgeGaussians) // Call the constructor of base class
{
    // Four corner nodes
    nodeCoord_[0] << -1, -1;
    nodeCoord_[1] << 1, -1;
    nodeCoord_[2] << 1, 1;
    nodeCoord_[3] << -1, 1;
}

ShapeQ4::~ShapeQ4()
{
}

MatrixXd ShapeQ4::functionMat(const Vector2d & point) const
{
    MatrixXd result = MatrixXd::Zero(2, 2 * numNodes_); // 2-D, 4-Node
    for (int i = 0; i < numNodes_; i++) {
        // Ni = (1+xi_i*xi)(1+eta_i*eta)/4
        double Ni = (1 + nodeCoord_[i](0) * point(0)) * (1 + nodeCoord_[0](1) * point(1)) / 4;
        result(0, 2 * i) = Ni;
        result(1, 2 * i + 1) = Ni;
    }
    return result;
}

MatrixXd ShapeQ4::functionDeriv(const Vector2d & point) const
{
    MatrixXd result(numNodes_, 2);
    for (int i = 0; i < numNodes_; i++) {
      // d(Ni)/d(xi)=xi_i*(1+eta_i*eta)/4
      result(i, 0) = nodeCoord_[i](0) * (1 + nodeCoord_[i](1) * point(1)) / 4;
      // d(Ni)/d(eta)=eta_i*(1+xi_i*xi)/4
      result(i, 1) = nodeCoord_[i](1) * (1 + nodeCoord_[i](0) * point(0)) / 4;
    }
    return result;
}
