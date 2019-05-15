/**
 * @file ViscoElastic.cpp
 * Implementation of ViscoElastic class.
 *
 * @author Haohang Huang111
 * @date May 19, 2018
 */

#include "ViscoElastic.h"

ViscoElastic::ViscoElastic(const bool & anisotropy, const bool & nonlinearity, const bool & noTension, const std::vector<double> & properties)
  : Material(anisotropy, nonlinearity, noTension)
{
    int i = 0;

    // Isotropic: Modulus, Poisson's ratio, body force (r,z), thermal coefficient, temperature change
    double K0 = properties[i++];
    double G0 = properties[i++];
    double t = properties[i++];

    double e = 2.718281828;

    // All the following are the temporary values.
    // G_array and K_array are supposed to be the same length.
    double G_array[] = {0.362,0.363,0.177,0.074,0.017,0.006};
    std::vector<double> G (G_array, G_array + sizeof(G_array) / sizeof(double) );
    double K_array[] = {0.362,0.363,0.177,0.074,0.017,0.006};
    std::vector<double> K (K_array, K_array + sizeof(K_array) / sizeof(double) );
    double tao_array[] = {4.09*0.000001,2.56*0.0001,7.71*0.001,2.1*0.1,3.88,6.53*10};
    std::vector<double> tao (tao_array, tao_array + sizeof(tao_array) / sizeof(double) );

    double Gt = G0;
    double Kt = K0;
    for (int j = 0; j < G.size(); j++){
        Gt -= G0*G[j]*(1 - std::pow(e, -t/tao[j]));
        Kt -= K0*K[j]*(1 - std::pow(e, -t/tao[j]));
    }
    M_ = 9*Kt*Gt/(3*Kt + Gt);
    v_ = (3*Kt - 2*Gt)/(6*Kt + 2*Gt);
    // Compute the stress-strain constitutive matrix
    E_ << 1 - v_, v_, v_, 0,
          v_, 1 - v_, v_, 0,
          v_, v_, 1 - v_, 0,
          0,  0,   0,  (1-2*v_)/2;
    E_ = E_ * M_ / (1+v_) /(1-2*v_);

    // Assign body force (unit weight)
    bodyForce_ << properties[i], properties[i+1];
    i += 2;

    // Assign thermal parameters
    double alpha  = properties[i++];
    double deltaT = properties[i++]; /** The temperature change (assume same across the entire element) */
    double strain = alpha * deltaT;
    thermalStrain_ << strain, strain, strain, 0;
}

ViscoElastic::~ViscoElastic()
{
}
