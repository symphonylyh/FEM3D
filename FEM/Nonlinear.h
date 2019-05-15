/**
 * @file Nonlinear.h
 * Derived class from Analysis for nonlinear elastic problems.
 *
 * @author Haohang Huang
 * @date May 19, 2018
 * @note Add no tension analysis on June 1, 2018.
 */

// Confusion:
// 1. principal stress compute in GT-PAVE
// 2. is my understanding of nonlinear analysis right?
// 3. the meaning of 3 Uzan parameters?
// 4. why the nonlinear converge to the same displacment result with linear elastic, but with very different (much smaller) modulus values?
// 5. the displacment values are related to the guess value?
#ifndef Nonlinear_h
#define Nonlinear_h

#include "Analysis.h"

/* Derived class for solving nonlinear elastic problems.
 */
class Nonlinear : public Analysis
{
  public:
    /* See the documentation of base class Analysis.
     */
    Nonlinear(Mesh & meshInfo);  // ctor cannot be inherited, should explicitly call base class's ctor in derived class's ctor
    ~Nonlinear();
    void solve();

    /**
     * Compute principal stresses at Gaussian points, and update the modulus and E matrix for next iterations.
     *
     * @return A boolean value incidating the convergence status at this iteration.
     */
    bool nonlinearIteration(double damping);

    /**
     * Compute unbalanced tension stresses at Gaussian points, and update the nodal force vector for next iterations.
     *
     * @return A boolean value incidating the convergence status at this iteration.
     */
    bool noTensionIteration();

    /**
     * Helper function for the convertion to principal stresses from cylindrical coordinates.
     *
     * @param stress Stresses in cylindrical coordinates, sigma_r, sigma_theta, sigma_z, tau_rz
     * @return The principal stresses in sigma3, sigma2, sigma1 order.
     */
    VectorXd principalStress(const VectorXd & stress) const;

  private:
    int gravityIncrementNum; /* No. of body load (gravity & temperature & residual) increments */
    int loadIncrementNum; /* No. of traffic load (point & edge) increments */
    double gravityDamping; /* Damping ratio lambda for body force incremental loading */
    double loadDamping; /* Damping ratio lambda for traffic incremental loading */

};

#endif /* Nonlinear_h */
