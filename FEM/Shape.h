/**
 * @file Shape.h
 * Abstract base class for various isoparametric element shapes.
 *
 * @author Haohang Huang
 * @date Feburary 10, 2018
 * @note Efficiency optimized by pass/return-by-ref of Gaussian points on March 25,
 * 2018
 * @note Efficiency optimized by pre-cached the shape functions evaluated at all
 * Gaussian integraton points on Apr 22, 2018.
 * @todo T3, Q4, Q6 shape functions to be implemented.
 */

#ifndef Shape_h
#define Shape_h

#include "Eigen/Eigen" // #include "Matrix/Eigen" if change name
#include <vector>

using namespace Eigen;

/* Abstract base Shape class with pre-cached shape functions evaluated at all
 * Gaussian integration points.
 *
 * @note For the polymorphism of a combination of different shapes (e.g., T3,
 * Q6, Q8), we create a base class that allows different implementation of specific
 * shapes using C++ inheritance features. The instance of this class is stored
 * as a static member variable in the associated Element class. For all the
 * instances of the associated Element class, only one single instance of Shape
 * class is allocated and can be accessed by pointer. Additionally, the shape
 * functions at Gaussian points are pre-cached and can be accessed in O(1)
 * running time. The reason of pre-caching is that our element and shape are
 * isoparametric, namely no matter how the global coordinates of nodes change,
 * the Gaussian points are always the same in local coordinates. Pre-caching
 * enables efficient queries of shape functions without repeatedly compute them.
 */
class Shape
{
    public:
        /**
         * Constructor for Shape class.
         *
         * @param nodes Number of nodes of this shape of element.
         * @param gaussians Number of Gaussian integration points of this
         * shape of element.
         * @param edges Number of edges of this shape of element.
         * @param edgeNodes Number of nodes each edge consists of.
         * @param edgeGaussians Number of Gaussian intergration points at each
         * edge.
         *
         * @note The information about edges and its integration are specified at
         * initialization because no assumption or generic forms can be made.
         * It should leave the decision for user to define and select different
         * types of integration scheme (normal integration, reduced integration,
         * etc.). Additionally, no assumption can be made because the shape could
         * be triangular, quadrilateral, or any arbitrary form so that no
         * predetermined Node number-Edge number-Gaussian point number relation
         * exists.
         */
        Shape(const int & nodes, const int & gaussians, const int & edges, const int & edgeNodes, const int & edgeGaussians);

        /**
         * Destructor for Shape.
         */
        virtual ~Shape() { }

        /**
         * Compute the shape Function N in a column vector at a given point.
         *
         * If n is the number of nodes of this element type, then this function will
         * returns a n-by-1 column vector:
         *
         * [ N1(x), N2(x), ..., Nn(x) ]
         *
         * where Ni(x) is the shape function with respect to ith node at point x.
         *
         * @param point A 2D point x(r,z) at which the shape function to be evaluated. Pass-by-ref.
         * @return The shape function N as a column vector at given point.
         *
         * @note Use a const qualifier for reference variable is the recommended programming style.
         * Shape function depends on given point, so we can only return-by-value, no
         * more space for optimization.
         */
        virtual VectorXd functionVec(const Vector2d & point) const = 0;

        /**
         * Stack and interleave the shape Function N above into a matrix for the
         * convenience of stiffness matrix computation.
         *
         * If n is the number of nodes of this element type, then this function will
         * returns a 2-by-2n matrix:
         *
         * [ N1(x)  0 | N2(x)  0 | ... | Nn(x)  0 ]
         * [ 0  N1(x) | 0  N2(x) | ... | 0  Nn(x) ]
         *
         * this is because our displacement vector is stored in an interleaving form
         * as [r1, z1, r2, z2, ..., rn, zn].
         *
         * @param point A 2D point x(r,z) at which the shape function to be evaluated. Pass-by-ref.
         * @return The stacked and interleaved shape function N as a matrix at given point.
         */
        virtual MatrixXd functionMat(const Vector2d & point) const = 0;

        /**
         * Compute the local derivatives dN/d(xi) and dN/d(eta) of the shape function
         * at a given point in isoparametric coordinates.
         *
         * If n is the number of nodes of this element type, then this function will
         * returns a 2-by-n matrix:
         *
         * [ dN1(x)/d(xi)   dN2(x)/d(xi)  ... dNn(x)/d(xi)  ]
         * [ dN1(x)/d(eta)  dN2(x)/d(eta) ... dNn(x)/d(eta) ]
         *
         * @param point A 2D point x(r,z) at which the shape function to be evaluated. Pass-by-ref.
         * @return The stacked and interleaved shape function derivatives dN as a matrix at given point.
         */
        virtual MatrixXd functionDeriv(const Vector2d & point) const = 0;

        /**
         * Compute the edge shape Function N in a column vector at a given point.
         *
         * @param point A 1D point x at which the edge shape function to be evaluated. Pass-by-ref.
         * @return The edge shape function N as a column vector at given point.
         */
        virtual VectorXd edgeFunctionVec(const double & point) const = 0;

        /**
         * Stack and interleave the edge shape Function N above into a matrix.
         *
         * @param point A 1D point x at which the edge shape function to be evaluated. Pass-by-ref.
         * @return The stacked and interleaved edge shape function N as a matrix at given point.
         */
        virtual MatrixXd edgeFunctionMat(const double & point) const = 0;

        /**
         * Compute the local derivatives dN/dx of the edge shape function
         * at a given point in isoparametric coordinates.
         *
         * @param point A 1D point x at which the edge shape function to be evaluated. Pass-by-ref.
         * @return The edge shape function derivatives dN as a column vector at given point.
         */
        virtual VectorXd edgeFunctionDeriv(const double & point) const = 0;

        /**
         * Get the gaussian integration points to be used in load distribution and
         * stress/strain computation.
         *
         * @return An array of 2D points. Since the whole Shape object is
         * dynamically allocated on heap, we can return-by-ref.
         */
        const std::vector<Vector2d> & gaussianPt() const;

        /**
         * Get the ith gaussian integration points.
         *
         * @param i The index of Gaussian points to be queried.
         * @return The Gaussian point.
         */
        const Vector2d & gaussianPt(const int & i) const;

        /**
         * Get the gaussian integration weights to be used in load distribution and
         * stress/strain computation.
         *
         * @return An array of weight values. Since the whole Shape object is
         * dynamically allocated on heap, we can return-by-ref.
         */
        const std::vector<double> & gaussianWt() const;

        /**
         * Get the ith gaussian weights.
         *
         * @param i The index of Gaussian points to be queried.
         * @return The Gaussian weight.
         */
        const double & gaussianWt(const int & i) const;

        /**
         * Get pre-cached shape function (vector form) at certin Gaussian point.
         * @param i The index of Gaussian points to be queried.
         * @return The shape function N as a column vector at Gaussian point i.
         *
         * @note Note that there are two sets of functions that have same name
         * but different function signatures (i.e., have different types of
         * parameters). The above one is used to compute shape function at any
         * arbitrary point, this function is to get the pre-cached shape function
         * at Gaussian points.
         */
        const VectorXd & functionVec(const int & i) const;

        /**
         * Get pre-cached shape function (matrix form) at certin Gaussian point.
         * @param i The index of Gaussian points to be queried.
         * @return The shape function N as a matrix (interleaved) at Gaussian point i.
         */
        const MatrixXd & functionMat(const int & i) const;

        /**
         * Get pre-cached shape function derivative matrix at certin Gaussian point.
         * @param i The index of Gaussian points to be queried.
         * @return The derivatives of shape function dN as a matrix (not interleaved) at Gaussian point i.
         */
        const MatrixXd & functionDeriv(const int & i) const;

        /**
         * Get pre-cached edge shape function (vector form) at certin Gaussian point.
         * @param i The index of Gaussian points to be queried.
         * @return The edge shape function N as a column vector at Gaussian point i.
         */
        const VectorXd & edgeFunctionVec(const int & i) const;

        /**
         * Get pre-cached edge shape function (matrix form) at certin Gaussian point.
         * @param i The index of Gaussian points to be queried.
         * @return The edge shape function N as a matrix (interleaved) at Gaussian point i.
         */
        const MatrixXd & edgeFunctionMat(const int & i) const;

        /**
         * Get pre-cached edge shape function derivative matrix at certin Gaussian point.
         * @param i The index of Gaussian points to be queried.
         * @return The derivatives of shape function dN as a column vector at Gaussian point i.
         */
        const VectorXd & edgeFunctionDeriv(const int & i) const;

        /**
         * Get the edge gaussian integration points to be used in load distribution and
         * stress/strain computation.
         *
         * @return An array of 1D points.
         */
        const std::vector<double> & edgeGaussianPt() const;

        /**
         * Get the edge gaussian integration weights to be used in load distribution and
         * stress/strain computation.
         *
         * @return An array of weight values.
         */
        const std::vector<double> & edgeGaussianWt() const;

        /**
         * Get the node index of ith edge.
         *
         * @param i The zero-based edge index being queried.
         * @return An array of node indices belonging to ith edge.
         */
        const std::vector<int> & edge(const int & i) const;

    protected:

        /**
         * Private helper function for Pre-caching the shape functions at Gaussian integration points.
         */
        void _cacheShape();

        // Example dimensions of the member variables are given for a Q8 element as follows

        /** The number of nodes of this element type */
        int numNodes_; // 8

        /** The local isoparametric coordinates (xi,eta) of the element nodes */
        std::vector<Vector2d> nodeCoord_; // length 8

        /** The number of Gaussian integration points of the element */
        int numGaussianPts_; // 9

        /** An array of the coordinates of 2D Gaussian integration points */
        std::vector<Vector2d> gaussianPt_; // length 9

        /** An array of the weights of Gaussian integration points */
        std::vector<double> gaussianWt_; // length 9

        /** An array of the shape functions (vector form) at all Gaussian points */
        std::vector<VectorXd> shapeVec_; // length 9, each 8x1 vector

        /** An array of the shape functions (stacked matrix form) at all Gaussian points */
        std::vector<MatrixXd> shapeMat_; // length 9, each 2x16 matrix (interleaved)

        /** An array of the derivatives of shape function (stacked matrix form) at all Gaussian points */
        std::vector<MatrixXd> shapeDeriv_; // length 9, each 2x8 matrix (not interleaved)

        /** The number of edges of this element type */
        int numEdges_; // 4

        /** The number of nodes an edge has */
        int numEdgeNodes_; // 3

        /** An array of the edges with its node indices */
        std::vector<std::vector<int> > edgeList_; // length 4, each length 3

        /** The number of Gaussian integration points of an edge */
        int numEdgeGaussianPts_; // 3

        /** An array of the coordinates of the 1D Gaussian points of an edge */
        std::vector<double> edgeGaussianPt_; // length 3

        /** An array of the weights of Gaussian points of an edge */
        std::vector<double> edgeGaussianWt_; // length 3

        /** An array of the edge shape functions (vector form) at Gaussian points */
        std::vector<VectorXd> edgeShapeVec_; // length 3, each 3x1 vector

        /** An array of the edge shape functions (stacked matrix form) at Gaussian points */
        std::vector<MatrixXd> edgeShapeMat_; // length 3, each 2x6 matrix (interleaved)

        /** An array of the derivatives of edge shape function (vector form) at Gaussian points */
        std::vector<VectorXd> edgeShapeDeriv_; // length 3, each 3x1 vector

};

#endif /* Shape_h */
