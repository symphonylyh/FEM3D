/**
 * @file NonlinearElastic.cpp
 * Implementation of NonlinearElastic class.
 *
 * @author Haohang Huang
 * @date May 19, 2018
 */

#include "NonlinearElastic.h"
#include <cmath>
#include <iostream>
NonlinearElastic::NonlinearElastic(const bool & anisotropy, const bool & nonlinearity, const bool & noTension, const std::vector<double> & properties, const int & model, const std::vector<double> & parameters)
  : Material(anisotropy, nonlinearity, noTension), modelNo(model), coeff(parameters)
{
    // Just copy from LinearElastic constructor
    int i = 0;
    if (!anisotropy) { // Isotropic: Modulus, Poisson's ratio, body force (r,z), thermal coefficient, temperature change
        M_ = properties[i++];
        v_ = properties[i++];

        // Compute the stress-strain constitutive matrix
        E_ << 1 - v_, v_, v_, 0,
              v_,   1-v_, v_, 0,
              v_,   v_,  1-v_, 0,
              0,  0,    0,  (1-2*v_)/2;
        E_ = E_ * M_ / (1+v_) /(1-2*v_);
    }
    else { // Cross-anisotropic: Modulus R, Modulus Z, Poisson's ratio R, Poisson's ratio Z, Shear Modulus G, body force (r,z), thermal coefficient, temperature change
        Mr_ = properties[i++];
        Mz_ = properties[i++];
        vr_ = properties[i++];
        vz_ = properties[i++];
        G_ = properties[i++];

        // Helper coefficients for filling E matrix
        double n = Mr_ / Mz_;
        double m = G_ / Mz_;
        double A = Mz_ / (1 + vr_) / (1 - vr_ - 2 * n * vz_ * vz_);

        // Compute the stress-strain constitutive matrix
        E_ << n * (1 - n * vz_ * vz_), n * (vr_ + n * vz_ * vz_), n * vz_ * (1 + vr_), 0,
              n * (vr_ + n * vz_ * vz_), n * (1 - n * vz_ * vz_), n * vz_ * (1 + vr_), 0,
              n * vz_ * (1 + vr_), n * vz_ * (1 + vr_), 1 - vr_ * vr_, 0,
              0, 0, 0, m * (1 + vr_) * (1 - vr_ - 2 * n * vz_ * vz_);
        E_ = E_ * A;
    }

    // Assign body force (unit weight)
    bodyForce_ << properties[i], properties[i+1];
    i += 2;

    // Assign thermal parameters
    double alpha  = properties[i++];
    double deltaT = properties[i++]; /** The temperature change (assume same across the entire element) */
    double strain = alpha * deltaT;
    thermalStrain_ << strain, strain, strain, 0;
}

NonlinearElastic::~NonlinearElastic()
{
}

VectorXd NonlinearElastic::stressDependentModulus(const VectorXd & stress) const
{
    // stress(0)-sigma3; stress(1)-sigma2; stress(2)-sigma1
    // Bulk stress: theta = sigma1 + sigma2 + sigma3
    // Deviatoric stress: sigma_d = sigma1 - sigma3
    // Octahedral shear stress (usually meaningful in 3D): tau_oct = sqrt((1-2)^2 + (2-3)^2 + (1-3)^2) / 3, for 2D := sqrt(2)/3 * sigma_d
    double bulk = std::abs( stress(2) + stress(1) + stress(0) );
    double deviator = std::abs( stress(2) - stress(0) );
    double octahedral = std::sqrt(std::pow(stress(2) - stress(1), 2) + std::pow(stress(1) - stress(0), 2) + std::pow(stress(2) - stress(0), 2)) / 3;
    double atm = 14.696; // atmospheric pressure, 101.325 kPa or 14.7 psi

    // Calculated new stress-dependent resilient modulus from different models
    double Mr = 0, Mz = 0, G = 0; // horizontal modulus, vertical modulus, shear modulus
    switch (modelNo) {
        case 1 :
            // Hicks and Monismith (1971) K-theta model
            // M = k1 * theta^k2
            if (!anisotropy) {
                Mz = coeff[0] * std::pow(bulk, coeff[1]);
            } else {
                Mr = coeff[0] * std::pow(bulk, coeff[1]); // although K-theta only has 2 parameters, we fill 0 to complement the triplet
                Mz = coeff[3] * std::pow(bulk, coeff[4]);
                G = coeff[6] * std::pow(bulk, coeff[7]);
            }
            break;
        case 2 :
            // Uzan (1985) model:
            // M = k1 * theta^k2 * sigma_d^k3
            if (!anisotropy) {
                Mz = coeff[0] * std::pow(bulk, coeff[1]) * std::pow(deviator, coeff[2]);
            } else {
                Mr = coeff[0] * std::pow(bulk, coeff[1]) * std::pow(deviator, coeff[2]);
                Mz = coeff[3] * std::pow(bulk, coeff[4]) * std::pow(deviator, coeff[5]);
                G = coeff[6] * std::pow(bulk, coeff[7]) * std::pow(deviator, coeff[8]);
            }
            break;
        case 3 :
            // Pezo (1993) UT-Austin model
            // M = k1 * sigma_d^k2 * sigma_3^k3
            if (!anisotropy) {
                Mz = coeff[0] * std::pow(deviator, coeff[1]) * std::pow(std::abs(stress(0)), coeff[2]);
            } else {
                Mr = coeff[0] * std::pow(deviator, coeff[1]) * std::pow(std::abs(stress(0)), coeff[2]);
                Mz = coeff[3] * std::pow(deviator, coeff[4]) * std::pow(std::abs(stress(0)), coeff[5]);
                G = coeff[6] * std::pow(deviator, coeff[7]) * std::pow(std::abs(stress(0)), coeff[8]);
            }
            break;
        case 4 :
            // NCHRP 1-37A (2004) MEPDG model
            // M = k1 * p_a * (theta/p_a)^k2 * (tau_oct/p_a + 1)^k3
            if (!anisotropy) {
                Mz = coeff[0] * atm * std::pow(bulk/atm, coeff[1]) * std::pow(octahedral/atm + 1, coeff[2]);
            } else {
                Mr = coeff[0] * atm * std::pow(bulk/atm, coeff[1]) * std::pow(octahedral/atm + 1, coeff[2]);
                Mz = coeff[3] * atm * std::pow(bulk/atm, coeff[4]) * std::pow(octahedral/atm + 1, coeff[5]);
                G = coeff[6] * atm * std::pow(bulk/atm, coeff[7]) * std::pow(octahedral/atm + 1, coeff[8]);
            }
            break;
    }
    VectorXd result(3);
    result << Mr, Mz, G;
    return result;
}

MatrixXd NonlinearElastic::EMatrix(const VectorXd & modulus) const
{
    MatrixXd E(4,4);
    if (!anisotropy) {
        E << 1 - v_, v_, v_, 0,
              v_,   1-v_, v_, 0,
              v_,   v_,  1-v_, 0,
              0,  0,    0,  (1-2*v_)/2;
        E = E * modulus(0) / (1+v_) /(1-2*v_);
    } else {
        double n = modulus(0) / modulus(1); // Mr / Mz
        double m = modulus(2) / modulus(1); // G / Mz
        double A = modulus(1) / (1 + vr_) / (1 - vr_ - 2 * n * vz_ * vz_);
        E <<  n * (1 - n * vz_ * vz_), n * (vr_ + n * vz_ * vz_), n * vz_ * (1 + vr_), 0,
              n * (vr_ + n * vz_ * vz_), n * (1 - n * vz_ * vz_), n * vz_ * (1 + vr_), 0,
              n * vz_ * (1 + vr_), n * vz_ * (1 + vr_), 1 - vr_ * vr_, 0,
              0, 0, 0, m * (1 + vr_) * (1 - vr_ - 2 * n * vz_ * vz_);
        E = E * A;
    }
    return E;
}
