/**
 * @file Analysis.h
 * Abstract Analysis class for various types of analysis.
 *
 * @author Haohang Huang
 * @date Feburary 26, 2018
 * @note Efficiency optimized by local stiffness matrix return-by-ref on March 27, 2018.
 * @note Efficiency optimized by the abstraction of all element- and shape-related operations
 * into generic forms on Apr 22, 2018.
 * @todo Query point method.
 */

#ifndef Analysis_h
#define Analysis_h

#include "Mesh.h"

/* Abstract base Analysis class with shared public methods and pure virtual methods.
 *
 * Key concepts of during the design of this class:
 *
 * 1. One of the benefits of C++ inheritance feature is the code can be kept
 * in the most concise way. For example, we defined a Analysis base class here,
 * and we will have derived class such as linear elastic class, nonlinear elastic,
 * viscoelastic, etc. Furthermore, we can also have classes derived from the
 * already derived class, like linear isotropic and linear anisotropic class from
 * the linear elastic class. The subroutines/functions/procedures that are
 * shared by the derived classes can be abstracted and defined in the base class,
 * while derived classes can still define other subroutines/functions/procedures
 * as their own features. In this way, our program is elegantly organized and
 * extendable, with the minimum redundancy in codes.
 *
 * 2. The information such as mesh, global stiffness, force, etc should be used
 * in any type of the problems, so they have been made public (specifically,
 * protected) in Analysis class. By inheritance we have all those info as public
 * in the derived classes. It is noted that we should avoid name confliction
 * when we define the own variables of the derived classes.
 * In addition, from a memory allocation view, the base class's memory space is
 * allocated inside the memory space of this derived class. But these space are
 * just copies from the base class, i.e., if we have a second derived class,
 * these space are independent and uniquely of their own.
 */
class Analysis
{
    public:
        /**
         * Default constructor for Analysis.
         */
        // Analysis();

        /**
         * Custom constructor to get the mesh info of input file for analysis.
         *
         * @param mesh The read-in mesh object.
         */
        // Analysis(std::string const & fileName);
        Analysis(Mesh & meshInfo);

        /**
         * Virtual destructor for abstract class.
         */
        virtual ~Analysis();

        /**
         * Assemble the global stiffness matrix from the local stiffness matrix
         * of each element, meanwhile assemble the force vector with modifications
         * based on the applied boundary conditons.
         */
        void assembleStiffness();

        /**
         * Apply point load and edge load at each node in the global force vector.
         * The body force and temperature load should be applied element-wise
         * during the stiffness matrix assembly steps.
         */
        void applyForce();

        /**
         * Designate boundary condition at nodes (not used).
         */
        // void boundaryCondition();

        /**
         * Solve the problem using different approaches.
         */
        virtual void solve() = 0;

        /**
         * Compute nodal strain and stress from the extrapolation of results at
         * Gaussian integration points and average the value over the mesh.
         */
        void computeStrainAndStress();

        /**
         * Average nodal strain and stress at each node.
         */
        void averageStrainAndStress();

        // @TODO Query queryPoint(const double & x, const double & y) const; // create a Query object of disp, strain, stress and return it. Query needs interpolation based on shape function, and need to find the nearest neighbor or locate the element it's in. Should also check if post-processing software does this for us

        /**
         * Print the nodal displacement.
         */
        void printDisp() const;

        /**
         * Print the nodal strain.
         */
        void printStrain() const;

        /**
         * Print the nodal stress.
         */
        void printStress() const;

        /**
         * Write the outputs to file.
         */
        void writeToFile(std::string const & fileName) const;

        /**
         * Write the outputs into a .vtk file for post-processing.
         */
        void writeToVTK(std::string const & fileName) const;

    protected: // "protected" is a good choice. Only visible to the derived class

        /** The mesh information of the problem */
        Mesh & mesh;

        /** The global stiffness matrix as a 2n-by-2n sparse matrix */
        SparseMatrix<double> globalStiffness;

        /** The nodal displacement 2n-by-1 vector */
        VectorXd nodalDisp;

        /** The nodal force 2n-by-1 vector */
        VectorXd nodalForce;

        /** The nodal strain n-by-4 matrix */
        MatrixXd nodalStrain;

        /** The nodal strain n-by-4 matrix */
        MatrixXd nodalStress;
};

#endif /* Analysis_h */
