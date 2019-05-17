/**
 * @file ShapeB20.h
 * Derived class from Shape for the shape information of isoparametric B20 element.
 *
 * @author Haohang Huang
 * @date Feburary 13, 2018
 * @note Efficiency optimized by pass/return-by-ref of Gaussian points on March 25,
 * 2018.
 * @note Efficiency optimized by pre-cached the shape functions evaluated at all
 * Gaussian integraton points on Apr 22, 2018.
 */

#ifndef ShapeB20_h
#define ShapeB20_h

#include "Shape.h"

/* Derived class for the shape of isoparametric B20 element.
 * The sketch and index of the B20 element is:
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
class ShapeB20 : public Shape
{
    public:
        /* See the documentation of base class Shape. */
        ShapeB20(const int & nodes, const int & gaussians, const int & edges, const int & edgeNodes, const int & edgeGaussians);
        ~ShapeB20();

        VectorXd functionVec(const Vector2d & point) const;
        MatrixXd functionMat(const Vector2d & point) const;
        MatrixXd functionDeriv(const Vector2d & point) const;

        VectorXd edgeFunctionVec(const double & point) const;
        MatrixXd edgeFunctionMat(const double & point) const;
        VectorXd edgeFunctionDeriv(const double & point) const;

};

#endif /* ShapeB20_h */
