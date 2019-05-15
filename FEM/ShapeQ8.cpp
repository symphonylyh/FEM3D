/**
 * @file ShapeQ8.cpp
 * Implementation of ShapeQ8 class.
 *
 * @author Haohang Huang
 * @date Feburary 13, 2018
 */

#include "ShapeQ8.h"
#include <cmath>

ShapeQ8::ShapeQ8(const int & nodes, const int & gaussians, const int & edges, const int & edgeNodes, const int & edgeGaussians) :
    Shape(nodes, gaussians, edges, edgeNodes, edgeGaussians) // Call the constructor of base class
{
    // Set up shape parameters
    // Local xi-eta coordinates of nodes
    // Corner nodes
    nodeCoord_[0] << -1, -1;
    nodeCoord_[1] << 1, -1;
    nodeCoord_[2] << 1, 1;
    nodeCoord_[3] << -1, 1;
    // Mid-side nodes
    nodeCoord_[4] << 0, -1;
    nodeCoord_[5] << 1, 0;
    nodeCoord_[6] << 0, 1;
    nodeCoord_[7] << -1, 0;

    // Local xi-eta coordinates of gaussian points
    double temp = std::sqrt(0.6);
    gaussianPt_[0] << -temp, -temp;
    gaussianPt_[1] << 0, -temp;
    gaussianPt_[2] << temp, -temp;
    gaussianPt_[3] << -temp, 0;
    gaussianPt_[4] << 0, 0;
    gaussianPt_[5] << temp, 0;
    gaussianPt_[6] << -temp, temp;
    gaussianPt_[7] << 0, temp;
    gaussianPt_[8] << temp, temp;
    // Edge Gaussian points
    edgeGaussianPt_[0] = -temp;
    edgeGaussianPt_[1] = 0;
    edgeGaussianPt_[2] = temp;

    // Gaussian weights
    double corner = 5.0 / 9.0; // do not use 5/9! that is 0 due to rounding down!
    double side = 8.0 / 9.0;
    // For the sake of better generailization in Element class, we directly compute
    // the 2D overlapped Gaussian weights at each Gaussian weights instead of
    // doing wi*wj later on. Thus, previously we need to:
    // for i = 1:3
    //     for j = 1:3
    //           result += X * w(i) * w(j);
    // Now we can just loop through all Gaussian points:
    // for g = 1:9
    //     result += X * W(g);
    // And the methods of Element class can be highly abstract.
    gaussianWt_[0] = corner * corner;
    gaussianWt_[2] = corner * corner;
    gaussianWt_[6] = corner * corner;
    gaussianWt_[8] = corner * corner;
    gaussianWt_[1] = side * corner;
    gaussianWt_[3] = side * corner;
    gaussianWt_[5] = side * corner;
    gaussianWt_[7] = side * corner;
    gaussianWt_[4] = side * side; // center
    // Edge Gaussian weights
    edgeGaussianWt_[0] = corner;
    edgeGaussianWt_[1] = side;
    edgeGaussianWt_[2] = corner;

    // Edge list
    std::vector<int> edge1{0,4,1};
    std::vector<int> edge2{1,5,2};
    std::vector<int> edge3{3,6,2};
    std::vector<int> edge4{0,7,3};
    edgeList_[0] = edge1;
    edgeList_[1] = edge2;
    edgeList_[2] = edge3;
    edgeList_[3] = edge4;

    // After setting up the above parameters, pre-cache the shape function
    _cacheShape();

}

ShapeQ8::~ShapeQ8()
{
}

VectorXd ShapeQ8::functionVec(const Vector2d & point) const
{ // 8x1 Vector
    VectorXd result(numNodes_);
    double Ni = 0;
    for (int i = 0; i < numNodes_; i++) {
        if (i < 4) { // 4 corner nodes
          // Ni = (1+xi_i*xi)(1+eta_i*eta)(xi_i*xi+eta_i*eta-1)/4
          Ni = (1 + nodeCoord_[i](0) * point(0)) * (1 + nodeCoord_[i](1) * point(1)) * (nodeCoord_[i](0) * point(0) + nodeCoord_[i](1) * point(1) - 1) / 4;
          // (Solved) @BUG here! Typo: wrote one nodeCoord_[i](0) to be nodeCoord_[0](0)
        }
        if (i == 4 || i == 6) { // xi = 0 mid-side nodes
          // Ni = (1-xi^2)(1+eta_i*eta)/2
          Ni = (1 - point(0) * point(0)) * (1 + nodeCoord_[i](1) * point(1)) / 2;
        }
        if (i== 5 || i == 7) { // eta = 0 mid-side nodes
          // Ni = (1+xi_i*xi)(1-eta^2)/2
          Ni = (1 + nodeCoord_[i](0) * point(0)) * (1 - point(1) * point(1)) / 2;
        }
        result(i) = Ni;
    }
    return result;
}

MatrixXd ShapeQ8::functionMat(const Vector2d & point) const
{ // 2x16 matrix
    MatrixXd result = MatrixXd::Zero(2, 2 * numNodes_); // 2-D, 8-Node
    double Ni = 0;
    for (int i = 0; i < numNodes_; i++) {
        if (i < 4) { // 4 corner nodes
          // Ni = (1+xi_i*xi)(1+eta_i*eta)(xi_i*xi+eta_i*eta-1)/4
          Ni = (1 + nodeCoord_[i](0) * point(0)) * (1 + nodeCoord_[i](1) * point(1)) * (nodeCoord_[i](0) * point(0) + nodeCoord_[i](1) * point(1) - 1) / 4;
        }
        if (i == 4 || i == 6) { // xi = 0 mid-side nodes
          // Ni = (1-xi^2)(1+eta_i*eta)/2
          Ni = (1 - point(0) * point(0)) * (1 + nodeCoord_[i](1) * point(1)) / 2;
        }
        if (i == 5 || i == 7) { // eta = 0 mid-side nodes
          // Ni = (1+xi_i*xi)(1-eta^2)/2
          Ni = (1 + nodeCoord_[i](0) * point(0)) * (1 - point(1) * point(1)) / 2;
        }
        result(0, 2 * i) = Ni;
        result(1, 2 * i + 1) = Ni;
    }
    return result;
}

MatrixXd ShapeQ8::functionDeriv(const Vector2d & point) const
{ // 2x8 matrix
    MatrixXd result(2, numNodes_);
    for (int i = 0; i < numNodes_; i++) {
      if (i < 4) { // 4 corner nodes
        // d(Ni)/d(xi)=xi_i*(1+eta_i*eta)*(2*xi_i*xi+eta_i*eta)/4
        result(0, i) = nodeCoord_[i](0) * (1 + nodeCoord_[i](1) * point(1)) * (2 * nodeCoord_[i](0) * point(0) + nodeCoord_[i](1) * point(1))/ 4;
        // d(Ni)/d(eta)=eta_i*(1+xi_i*xi)*(2*eta_i*eta+xi_i*xi)/4
        result(1, i) = nodeCoord_[i](1) * (1 + nodeCoord_[i](0) * point(0)) * (2 * nodeCoord_[i](1) * point(1) + nodeCoord_[i](0) * point(0))/ 4;
      }
      if (i == 4 || i == 6) { // xi = 0 mid-side nodes
        // d(Ni)/d(xi)= -xi*(1+eta_i*eta)
        result(0, i) = - point(0) * (1 + nodeCoord_[i](1) * point(1));
        // d(Ni)/d(eta)= (1-xi^2)*eta_i/2
        result(1, i) = (1 - point(0) * point(0)) * nodeCoord_[i](1) / 2;
      }
      if (i== 5 || i == 7) { // eta = 0 mid-side nodes
        // d(Ni)/d(xi)= (1-eta^2)*xi_i/2
        result(0, i) = (1 - point(1) * point(1)) * nodeCoord_[i](0) / 2;
        // d(Ni)/d(eta)= -eta*(1+xi_i*xi)
        result(1, i) = - point(1) * (1 + nodeCoord_[i](0) * point(0));
      }
    }
    return result;
}

VectorXd ShapeQ8::edgeFunctionVec(const double & point) const
{ // 3x1 vector
    VectorXd result(numEdgeNodes_);
    double Ni = 0;
    for (int i = 0; i < 3; i++) {
        if (i == 0) // left/bottom point
          Ni = - point * (1 - point) / 2; // Ni = -x(1-x)/2
        if (i == 2) // right/top point
          Ni = point * (1 + point) / 2; // Ni = x(1+x)/2
        if (i == 1) // mid-side point
          Ni = 1 - point * point; // Ni = (1-x^2)
        result(i) = Ni;
    }
    return result;
}

MatrixXd ShapeQ8::edgeFunctionMat(const double & point) const
{ // 2x6 matrix, 3 gaussian point and 3 node at each edge
    MatrixXd result = MatrixXd::Zero(2, 2 * numEdgeNodes_);
    double Ni = 0;
    for (int i = 0; i < 3; i++) {
        if (i == 0) // left/bottom point
          Ni = - point * (1 - point) / 2; // Ni = -x(1-x)/2
        if (i == 2) // right/top point
          Ni = point * (1 + point) / 2; // Ni = x(1+x)/2
        if (i == 1) // mid-side point
          Ni = 1 - point * point; // Ni = (1-x^2)
        result(0, 2 * i) = Ni;
        result(1, 2 * i + 1) = Ni;
    }
    return result;
}

VectorXd ShapeQ8::edgeFunctionDeriv(const double & point) const
{ // 3x1 vector
    VectorXd result(numEdgeNodes_);
    double Ni = 0;
    for (int i = 0; i < 3; i++) {
        if (i == 0) // left/bottom point
          Ni = (2 * point - 1) / 2; // dNi/dx = (2x-1)/2
        if (i == 2) // right/top point
          Ni = (2 * point + 1) / 2; // dNi/dx = (2x+1)/2
        if (i == 1) // mid-side point
          Ni = - 2 * point; // dNi/dx = -2x
        result(i) = Ni;
    }
    return result;
}
