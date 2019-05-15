// /**
//  * @file ShapeT6.h
//  * Derived class from Shape for the shape information of isoparametric T6 element.
//  *
//  * @author Haohang Huang
//  * @date Apr 22, 2018
//  */
//
// #ifndef ShapeT6_h
// #define ShapeT6_h
//
// #include "Shape.h"
//
// /* Derived class for the shape of isoparametric T6 element (a isosceles right triangle).
//  * The sketch and index of the T6 element is:
//  *
//  * 1 _
//  * |  \_
//  * 4    3_
//  * |      \_
//  * 2 -- 5 -- 0
//  *
//  * The sketch and index of the Gaussian integration points is:
//  *
//  * 5 _
//  * |  \_
//  * 3 -- 4 _
//  * |    |  \_
//  * 0 -- 1 -- 2
//  */
// class ShapeT6 : public Shape
// {
//     public:
//         /* See the documentation of base class Shape. */
//         ShapeT6(const int & nodes, const int & gaussians, const int & edges, const int & edgeNodes, const int & edgeGaussians);
//         ~ShapeT6();
//
//         VectorXd functionVec(const Vector2d & point) const;
//         MatrixXd functionMat(const Vector2d & point) const;
//         MatrixXd functionDeriv(const Vector2d & point) const;
//
//         VectorXd edgeFunctionVec(const double & point) const;
//         MatrixXd edgeFunctionMat(const double & point) const;
//         VectorXd edgeFunctionDeriv(const double & point) const;
//
// };
//
// #endif /* ShapeT6_h */
