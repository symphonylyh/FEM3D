/**
 * @file Linear.cpp
 * Implementation of Linear class.
 *
 * @author Haohang Huang
 * @date Feburary 26, 2018
 */

#include "Linear.h"
#include <chrono>
#include <iostream>

Linear::Linear(Mesh & meshInfo) : Analysis(meshInfo)
{
}

Linear::~Linear()
{
}

void Linear::solve()
{
    applyForce();
    assembleStiffness();
    // auto start = std::chrono::high_resolution_clock::now();
    //     assembleStiffnessAndForce();
    // auto finish = std::chrono::high_resolution_clock::now();
    // auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
    // std::cout << "Elapsed time assemble: " << elapsed.count() << " ms" << std::endl;

    // // All integrated into assembleStiffnessAndForce()
    // start = std::chrono::high_resolution_clock::now();
    //     applyForce();
    // finish = std::chrono::high_resolution_clock::now();
    // elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
    // std::cout << "Elapsed time apply force: " << elapsed.count() << " ms" << std::endl;
    //
    // start = std::chrono::high_resolution_clock::now();
    //     boundaryCondition();
    // finish = std::chrono::high_resolution_clock::now();
    // elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
    // std::cout << "Elapsed time apply boundary: " << elapsed.count() << " ms" << std::endl;

    // Option1: SimplicialLDLT <SparseMatrix<double> > solver;
    // Option2: ConjugateGradient <SparseMatrix<double> > solver;
    // SimplicialLDLT is direct solver: Recommended for very sparse and not too large problems (e.g., 2D Poisson eq.)
    // ConjugateGradient is iterative solver: Recommended for large symmetric problems (e.g., 3D Poisson eq.)
    // Ref: https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html
    // SimplicialLDLT <SparseMatrix<double> > solver;
    // ConjugateGradient <SparseMatrix<double> > solver;

    SimplicialLDLT <SparseMatrix<double> > solver;
    solver.compute(globalStiffness);
    nodalDisp = solver.solve(nodalForce);
    //
    // start = std::chrono::high_resolution_clock::now();
    //     // SparseLU <SparseMatrix<double> > solver;
    //     SimplicialLDLT <SparseMatrix<double> > solver;
    //     solver.compute(globalStiffness);
    // finish = std::chrono::high_resolution_clock::now();
    // elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
    // std::cout << "Elapsed time set sparse: " << elapsed.count() << " ms" << std::endl;
    //
    // start = std::chrono::high_resolution_clock::now();
    //     nodalDisp = solver.solve(nodalForce);
    // finish = std::chrono::high_resolution_clock::now();
    // elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
    // std::cout << "Elapsed solve sparse: " << elapsed.count() << " ms" << std::endl;

    // Compute strain and stress and accumulate at each node
    computeStrainAndStress();

    // Average strain and stress at each node
    averageStrainAndStress();
}
