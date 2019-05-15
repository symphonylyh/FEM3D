/**
 * @file BackAnalysis.cpp
 * Implementation of BackAnalysis class.
 *
 * @author Haohang Huang
 * @date September 16, 2018
 */

#include <iostream>
#include <cmath>
#include "BackAnalysis.h"
#include "Linear.h"

BackAnalysis::BackAnalysis(Mesh & meshInfo) : Analysis(meshInfo)
{
}

BackAnalysis::~BackAnalysis()
{
}

void BackAnalysis::solve()
{
    // Data of seven FWD sensors
    VectorXi sensorNode(7);
    VectorXd sensorDisp(7);
    sensorNode << 0, 12, 24, 36, 48, 60, 72; // user input
    sensorDisp << 0.01142, 0.00716,	0.00532,	0.00374,	0.00252,	0.00169,	0.00116; // user input
    double criteria = 0.005; // user input
    double precision = 0.5; // user input

    bool match = false;
    double compare = 0;
    double deviation = 0;
    // while (!match) {
    for (int x = 0; x < 10; x++) {
        // Adjust material properties for current iteration
        // What should be the feedback mechanism for iteration?
        // If computed displacement < sensor displacement, meaning the layer modulus are too high-->decrease; otherwise, increase

        // Material information is stored in Mesh class
        for (auto & m : mesh.materialList) {
            m->adjustModulus(compare * precision); // Note: I add this method for Material{.h,.cpp} and LinearElastic{.h,.cpp}. Don't forget.
        }
        // computed < measured, compare = -1-->layer modulus should decrease, so + (-1) * modulus. Consistent

        // Run analysis
        applyForce();
        assembleStiffness();

        SimplicialLDLT <SparseMatrix<double> > solver;
        solver.compute(globalStiffness);
        nodalDisp = solver.solve(nodalForce);

        // Check convergence
        match = true;
        compare = 0;
        for (int i = 0; i < static_cast<int>(sensorNode.size()); i++) {
            double error = - nodalDisp(2 * sensorNode(i) + 1) - sensorDisp(i);
            std::cout << "Sensor " << i << " Error = " << error << std::endl;
            deviation += (std::abs(error) / sensorDisp(i)) * 100;
            if ( std::abs(error) > criteria ) {
                // FWD only requires displacement value match. So strain and stress computation is unnecessary between iterations.
                match = false;
            }
            compare += error < 0 ? -1 : 1; // computed < measured, -1; computed > measured, 1
        }
        compare = compare < 0 ? -1 : 1; // 7 sensors, compare can never be tied as 0. If more sensors are at computed < measured, set compare to -1, and vice versa
        std::cout << "Modulus Adjust = " << compare << " (-1 decrease; 1 increase)" << std::endl;
        std::cout << "Converge = " << match << " (0 fail; 1 success)" <<std::endl;
        std::cout << "Average Error (%) = " << deviation / 7 << std::endl;
        std::cout << std::endl;
        deviation = 0;

    } // exit loop when vertical displacement at all sensor nodes are matching

    // Compute strain and stress and accumulate at each node
    computeStrainAndStress();

    // Average strain and stress at each node
    averageStrainAndStress();
}
