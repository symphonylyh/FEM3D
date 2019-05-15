/**
 * @file ShapeQ8.h
 * Derived class from Shape for the shape information of isoparametric Q8 element.
 *
 * @author Haohang Huang
 * @date Feburary 13, 2018
 * @note Efficiency optimized by pass/return-by-ref of Gaussian points on March 25,
 * 2018.
 * @note Efficiency optimized by pre-cached the shape functions evaluated at all
 * Gaussian integraton points on Apr 22, 2018.
 */

#ifndef ShapeQ8_h
#define ShapeQ8_h

#include "Shape.h"

/* Derived class for the shape of isoparametric Q8 element.
 * The sketch and index of the Q8 element is:
 *
 * 3 -- 6 -- 2
 * |         |
 * 7         5
 * |         |
 * 0 -- 4 -- 1
 *
 * The sketch and index of the Gaussian integration points is:
 *
 * 6 -- 7 -- 8
 * |    |    |
 * 3 -- 4 -- 5
 * |    |    |
 * 0 -- 1 -- 2
 */
class ShapeQ8 : public Shape
{
    public:
        /* See the documentation of base class Shape. */
        ShapeQ8(const int & nodes, const int & gaussians, const int & edges, const int & edgeNodes, const int & edgeGaussians);
        ~ShapeQ8();

        VectorXd functionVec(const Vector2d & point) const;
        MatrixXd functionMat(const Vector2d & point) const;
        MatrixXd functionDeriv(const Vector2d & point) const;

        VectorXd edgeFunctionVec(const double & point) const;
        MatrixXd edgeFunctionMat(const double & point) const;
        VectorXd edgeFunctionDeriv(const double & point) const;

};

#endif /* ShapeQ8_h */
