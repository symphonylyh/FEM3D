/**
 * @file ShapeQ4.h
 * Derived class from Shape for the shape information of isoparametric Q4 element.
 *
 * @author Haohang Huang
 * @date Feburary 10, 2018
 * @note Efficiency optimized by pass/return-by-ref of Gaussian points on March 25,
 * 2018
 * @note Efficiency optimized by pre-cached the shape functions evaluated at all
 * Gaussian integraton points on Apr 22, 2018.
 * @todo Gaussian points and functionVec not implemented yet.
 */

#ifndef ShapeQ4_h
#define ShapeQ4_h

#include "Shape.h"

/* Derived class for the shape of isoparametric Q4 element.
 * The sketch and index of the Q4 element is:
 *
 * 3 -- 2
 * |    |
 * 0 -- 1
 *
 * The sketch and index of the Gaussian integration points is:
 *
 * 3 -- 2
 * |    |
 * 0 -- 1
 */
class ShapeQ4 : public Shape
{
    public:
        /* See the documentation of base class Shape. */
        ShapeQ4(const int & nodes, const int & gaussians, const int & edges, const int & edgeNodes, const int & edgeGaussians);
        ~ShapeQ4();
        MatrixXd functionMat(const Vector2d & point) const;
        MatrixXd functionDeriv(const Vector2d & point) const;

};

#endif /* ShapeQ4_h */
