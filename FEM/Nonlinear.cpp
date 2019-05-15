/**
 * @file Nonlinear.cpp
 * Implementation of Nonlinear class.
 *
 * @author Haohang Huang
 * @date May 19, 2018
 */

#include "Nonlinear.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <functional>

Nonlinear::Nonlinear(Mesh & meshInfo) : Analysis(meshInfo)
{
    gravityIncrementNum = (mesh.iterations)[0];
    loadIncrementNum = (mesh.iterations)[1];
    gravityDamping = (mesh.iterations)[2];
    loadDamping = (mesh.iterations)[3];
}

Nonlinear::~Nonlinear()
{
}

void Nonlinear::solve()
{
// Several options that can affect the results:
// 1. Whether to do one last run after the convergence
// 2. Change damping ratio
// 3. Whether to do incremental loading, and how many increments
// 4. A strict (check every individual Gaussian) or loose (check only at center) convergence criteria, also to use different modulus for every Gaussian or for every element
// 5. Convergence criteria to compare with old modulus or new modulus
// no increment, all Gaussian but only check convergence at center, 0.3 damping, no final run --> 2.57e-4 result

bool incremental = true;
if (gravityIncrementNum == 0 && loadIncrementNum == 0) incremental = false; // input two zeros means I don't want to have incremental loading

if (incremental) {
    // -----------------------------------------------------------------------------
    // --------------- Start of Incremental Loading Scheme -------------------------
    // -----------------------------------------------------------------------------

    // Gravity, temperature, and residual stress increments
    // Idea: for each material, re-assign the body force, residual stress and
    // thermal strain incrementally. At the beginning we should calculate the increments
    // A good observation: with gravity load only, the stress is independent with the modulus,
    // so any arbitrary initial guess of the modulus won't affect the stress-dependent modulus.
    const std::vector<Material*> & materials = mesh.materialList;
    std::vector<Vector2d> gravityIncrement;
    std::vector<VectorXd> thermalIncrement;
    gravityIncrement.reserve(materials.size());
    thermalIncrement.reserve(materials.size());

    // Record the total gravity load
    for (auto & m : materials) {
        gravityIncrement.push_back(m->bodyForce() / gravityIncrementNum);
        thermalIncrement.push_back(m->thermalStrain() / gravityIncrementNum);
    }

    for (int ic = 1; ic <= gravityIncrementNum; ic++) {
        // Apply the load incrementally
        for (unsigned m = 0; m < materials.size(); m++) {
            materials[m]->setBodyForce(gravityIncrement[m] * ic);
            materials[m]->setThermalStrain(thermalIncrement[m] * ic);
            // std::cout << "Debug: " << gravityIncrement[m] * ic << std::endl;
        }

        // Achieve the modulus and tension convergence at each increment
        bool nonlinearConvergence = false;
        int count = 0;
        while (!nonlinearConvergence) { // convergence criteria
            // Assemble the K and F based on the mesh information (without applying any
            // load, this is for body force and temperature load incremental only)
            // @BUG(solved) Normally applyForce() will initialize the global force vector
            // and the assembleStiffness() function below will always do += for nodal force. But
            // in the body force increments, applyForce() is not called, therefore we
            // need to manually reset the nodalForce otherwise it will keeps accumulating.
            nodalForce = VectorXd::Zero(2 * mesh.nodeCount());
            assembleStiffness();

            // Solve K U = F
            SimplicialLDLT <SparseMatrix<double> > solver;
            solver.compute(globalStiffness);
            nodalDisp = solver.solve(nodalForce);

            // Traverse each element, compute stress at Gaussian points, and update the modulus for the next (i + 1) iteration (if current iteration is i)
            nonlinearConvergence = nonlinearIteration(gravityDamping);

            count++;
        }
        // For the exit iteration, the new converged modulus is updated, but the nodalDisp
        // is for the last iteration, so we should do one more solve to match the modulus & displacment
        nodalForce = VectorXd::Zero(2 * mesh.nodeCount());
        assembleStiffness();
        SimplicialLDLT <SparseMatrix<double> > solver;
        solver.compute(globalStiffness);
        nodalDisp = solver.solve(nodalForce);

        std::cout << "Body Force Increment No." << ic << ", Total iterations = " << count << std::endl;
        // std::cout << "Nodal Displacement: ";
        // std::cout << std::endl;
        // for (int i = 0; i < mesh.nodeCount(); i++) {
        //   std::cout << "Node " << i << " : " << nodalDisp(2 * i) << " " << nodalDisp(2 * i + 1) << std::endl;
        // }
        // std::cout << std::endl;
        std::cout << "-----------------------------------------" << std::endl;
    }
    // std::cout << "Material load applied! \n" << std::endl;
    // For triaxial case: Output the displacment information after the body load but before the surface load
    // std::cout << nodalDisp(2 * 28 + 1) << " " << nodalDisp(2 * 72 + 1) << std::endl;
    // OR
    // computeStrainAndStress();
    // averageStrainAndStress();
    // std::cout << mesh.nodeArray()[28]->getDisp()(1) << " " << mesh.nodeArray()[72]->getDisp()(1) << std::endl;

    // Should note the possible numerical error by dividing the increments, e.g.,
    // x / 5 * 5 might not be exactly the same number. Actually they are almost the same

    // Traffic load increments (point load and edge load)
    std::vector<double> pointLoadIncrement = mesh.loadValue;
    std::vector<std::vector<double> > edgeLoadIncrement = mesh.edgeLoadValue;
    std::transform(pointLoadIncrement.begin(), pointLoadIncrement.end(), pointLoadIncrement.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, 1.0 / loadIncrementNum)); // in-place change
    for (auto & e : edgeLoadIncrement)
        std::transform(e.begin(), e.end(), e.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, 1.0 / loadIncrementNum));

    for (int ic = 1; ic <= loadIncrementNum; ic++) {
        // Apply the load incrementally
        std::transform(pointLoadIncrement.begin(), pointLoadIncrement.end(), mesh.loadValue.begin(), std::bind(std::multiplies<int>(), std::placeholders::_1, ic));
        for (unsigned e = 0; e < edgeLoadIncrement.size(); e++)
            std::transform(edgeLoadIncrement[e].begin(), edgeLoadIncrement[e].end(), (mesh.edgeLoadValue)[e].begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, ic));

        // Achieve the modulus and tension convergence at each increment
        bool nonlinearConvergence = false;
        int count = 0;
        while (!nonlinearConvergence) { // convergence criteria
            // Assemble the K and F based on the mesh information (with traffic load applied)
            applyForce();
            assembleStiffness();

            // Solve K U = F
            SimplicialLDLT <SparseMatrix<double> > solver;
            solver.compute(globalStiffness);
            nodalDisp = solver.solve(nodalForce);

            // Traverse each element, compute stress at Gaussian points, and update the modulus for the next (i + 1) iteration (if current iteration is i)
            nonlinearConvergence = nonlinearIteration(loadDamping);

            count++;
        }
        // For the exit iteration, the new converged modulus is updated, but the nodalDisp
        // is for the last iteration, so we should do one more solve to match the modulus & displacment
        applyForce();
        assembleStiffness();
        SimplicialLDLT <SparseMatrix<double> > solver;
        solver.compute(globalStiffness);
        nodalDisp = solver.solve(nodalForce);

        std::cout << "Traffic Load Increment No." << ic << ", Total iterations = " << count << std::endl;
        std::cout << "-----------------------------------------" << std::endl;
    }
    // std::cout << "Traffic load applied! \n" << std::endl;

    // -----------------------------------------------------------------------------
    // ----------------- End of Incremental Loading Scheme -------------------------
    // -----------------------------------------------------------------------------
}
else { // experimental module for tension modification (loads are applied directly rather than incrementally for convenience)
    // -------------------------------------------------------------------------
    // --------------- Start of Nonlinear Iteration Scheme ---------------------
    // -------------------------------------------------------------------------
    bool nonlinearConvergence = false;
    int i = 0; // for debug print only
    while (!nonlinearConvergence) { // convergence criteria
    //for (int i = 0; i < 10; i++) {
        // std::cout << "Nonlinear Iteration No." << i++ << std::endl;
        // Assemble the K and F based on the mesh information. At 1st iteration, the initial guess modulus M0 will be used; later on at iteration i, the stress-dependent modulus updated from (i - 1) iteratiion will be used
        applyForce();
        assembleStiffness();

        // Solve K U = F
        SimplicialLDLT <SparseMatrix<double> > solver;
        solver.compute(globalStiffness);
        nodalDisp = solver.solve(nodalForce);

        // Traverse each element, compute stress at Gaussian points, and update the modulus for the next (i + 1) iteration (if current iteration is i)
        nonlinearConvergence = nonlinearIteration(0.3);
    }
    applyForce();
    assembleStiffness();
    SimplicialLDLT <SparseMatrix<double> > solver1;
    solver1.compute(globalStiffness);
    nodalDisp = solver1.solve(nodalForce);
    // After convergence is achieved at the last iteration, the solved displacment
    // is stored in the protected member of Analysis class -- nodalDisp. And
    // globalStiffness & nodalForce are also pre-cached. K, U, F are all knowns
    // and can be used in the following no tension iteration scheme.
    // -------------------------------------------------------------------------
    // -------------------- End of Nonlinear Scheme ----------------------------
    // -------------------------------------------------------------------------

    bool tensionScheme = false;
    if (tensionScheme) {
    // -------------------------------------------------------------------------
    // --------------- Start of No Tension Iteration Scheme --------------------
    // -------------------------------------------------------------------------
    bool tensionConvergence = false;
    SimplicialLDLT <SparseMatrix<double> > solver;
    solver.compute(globalStiffness);
    i = 0;
    while (!tensionConvergence) { // convergence criteria
    // for (int i = 0; i < 2; i++) { // for debug print only
        // std::cout << "Tension Iteration No." << i++ << std::endl;
        // Solve K U = F
        // Note 1: the above nonlinear iteration scheme is iteratively solving a
        // series of linear elastic cases, where K should be updated every time.
        // But in the no tension iteration scheme, K remains unchanged as the
        // last iteration in the nonlinear process. We only update the F vector.
        // Note 2: in Eigen, the solver.compute() is a pre-conditioning of matrix,
        // and we can just recycle the solver for current use. Therefore, the
        // solver is placed outside the while loop.
        nodalDisp = solver.solve(nodalForce);

        // Traverse each element, compute stress at Gaussian points, and update the modulus for the next (i + 1) iteration (if current iteration is i)
        tensionConvergence = noTensionIteration();
        std::cout << "Converge?: " << (tensionConvergence == true ? "Yes" : "No") << std::endl;
        // Thinking: we can add a check for divergence and stop the no-tension scheme, saying that the mesh/material is not suitable for no-tension modification
    }
    // -------------------------------------------------------------------------
    // ----------------  End of No Tension Iteration Scheme --------------------
    // -------------------------------------------------------------------------
    }
}
    // After both material nonlinearity and granular no-tension scheme converge,
    // compute the nodal strain and stress from the final displacment results
    computeStrainAndStress();
    averageStrainAndStress();

    // Output the average axial strain at the surface (for fastcell case only)
    int node_start = 0, node_end = 38 + 1; // Node 0 ~ 40 are the surface nodes
    int node_curr = node_start;
    double disp = 0;
    while (node_curr != node_end) {
        disp += mesh.nodeArray()[node_curr]->getDisp()(1);
        node_curr++;
    }
    double strain = - disp / (node_end - node_start) / (150.0/25.4);
    std::cout << strain << std::endl;
}

bool Nonlinear::nonlinearIteration(double damping)
{
    bool convergence = true;
    double sumError = 0;
    double sumModulus = 0;

    Element* curr;
    int numNodes; // number of nodes belong to the element
    int numGaussianPt; // number of Gaussian points of the element
    for (int i = 0; i < mesh.elementCount(); i++) {
        curr = mesh.elementArray()[i];
        Material* material = curr->material();
        if (material->nonlinearity) { // compute stress for nonlinear elastic element only, skip all linear elastic ones
            const VectorXi & nodeList = curr->getNodeList();
            numNodes = curr->getSize();
            numGaussianPt = (int)curr->shape()->gaussianPt().size();

            // Assemble the nodal displacement vector for this element
            VectorXd nodeDisp(2 * numNodes);
            for (int j = 0; j < numNodes; j++) {
                nodeDisp(2 * j) = nodalDisp(2 * nodeList(j));
                nodeDisp(2 * j + 1) = nodalDisp(2 * nodeList(j) + 1);
            }

            // Step 1: Compute stress at gaussian points based on cached M & E from last iteration
            // Step 2: Update new modulus based on the stress from step 1 and mix with old modulus via damping ratio
            // Step 3: Cache the modulus to be used in the next iteration
            // Note: Tutu's approach only use the center Gaussian point for the whole element, as follows
            // MatrixXd B = curr->BMatrix(curr->shape()->gaussianPt(4));
            // VectorXd strain = B * nodeDisp; // e = B * u
            // double modulus_old = (curr->modulusAtGaussPt)(4); // M_(i-1)
            // VectorXd modulus_old_vec << modulus_old;
            // VectorXd stress = material->EMatrix(modulus_old_vec) * (strain - curr->thermalStrain()); // sigma = E_(i-1) * (e - e0), note that the M and E are both from previous iteration
            // VectorXd modulus_vec = material->stressDependentModulus(principalStress(stress));
            // double modulus_new = modulus_vec(1); // M_i, 1 for vertical modulus
            // double modulus = (1 - damping) * modulus_old + damping * modulus_new; // true M_i after applying damping ratio

            // For isotropy case (or a simplified anisotropy case), iterate only on the single modulus (vertical modulus for anisotropy)
            if (!material->anisotropy) {
                for (int g = 0; g < numGaussianPt; g++) {
                    // More strict approach
                    MatrixXd B = curr->BMatrix(curr->shape()->gaussianPt(g));
                    VectorXd strain = B * nodeDisp; // e = B * u
                    double modulus_old = (curr->modulusAtGaussPt)(g); // M_(i-1)
                    VectorXd modulus_old_vec(1);
                    modulus_old_vec << modulus_old;
                    VectorXd stress = material->EMatrix(modulus_old_vec) * (strain - curr->thermalStrain()); // sigma = E_(i-1) * (e - e0), note that the M and E are both from previous iteration
                    // tension modification
                    // VectorXd principal = principalStress(stress);
                    // for (int x = 0; x < 3; x++) {
                    //     if (principal(x) > 0) principal(x) = 0;
                    // }
                    // VectorXd modulus_vec = material->stressDependentModulus(principal);
                    // double modulus_new = modulus_vec(1);
                    VectorXd modulus_vec = material->stressDependentModulus(principalStress(stress));
                    double modulus_new = modulus_vec(1); // M_i, 1 for vertical modulus
                    double modulus = (1 - damping) * modulus_old + damping * modulus_new; // true M_i after applying damping ratio

                    (curr->modulusAtGaussPt)(g) = modulus;

                    // Convergence criteria
                    // Criteria 1: modulus stabilize within 5% at all Gaussian points (less strict criteria only checks the center Gaussian point)
                    double error = std::abs(modulus - modulus_new);
                    if (g == 4 && error / modulus_old > 0.05) // tutu uses modulus_old, but I want to use modulus
                        convergence = false;
                        // should I just return false? No. Because you want to update all elements' modulus synchronously
                    // Criteraia 2: Accumulative modulus error within 0.2%
                    if (g == 4/*true*/) { // less strict convergence criteria
                        sumError += error * error;
                        sumModulus += modulus_old * modulus_old; // tutu uses modulus_old, but I want to use modulus
                    }
                    // For Debug Use
                    if (i == 1 && g == 4) { // the granular element at centerline
                        // std::cout << "nodelDisp: " << nodeDisp.transpose() << std::endl;
                        // std::cout << "Strain: " << strain.transpose() << std::endl;
                        // std::cout << "E: " << material->EMatrix(modulus_old) << std::endl;
                        // std::cout << "cylindrical stress: " << stress.transpose() << std::endl;
                        // std::cout << "principal stress: " << principalStress(stress).transpose() << std::endl;

                        // std::cout << "Old modulus: " << modulus_old << std::endl;
                        // std::cout << "New modulus: " << modulus_new << std::endl;
                        // std::cout << "True modulus: " << modulus << std::endl;
                    }
                }
            }
            // For anisotropy case, iterate on all 3 moduli (vertical, horizontal, shear modulus)
            else {
                for (int g = 0; g < numGaussianPt; g++) {
                    // More strict approach
                    MatrixXd B = curr->BMatrix(curr->shape()->gaussianPt(g));
                    VectorXd strain = B * nodeDisp; // e = B * u
                    VectorXd modulus_old = (curr->modulusAtGaussPt).row(g); // M_(i-1)
                    VectorXd stress = material->EMatrix(modulus_old) * (strain - curr->thermalStrain()); // sigma = E_(i-1) * (e - e0), note that the M and E are both from previous iteration
                    VectorXd modulus_new = material->stressDependentModulus(principalStress(stress)); // M_i
                    VectorXd modulus = (1 - damping) * modulus_old + damping * modulus_new; // true M_i after applying damping ratio

                    (curr->modulusAtGaussPt).row(g) = modulus;

                    // Convergence criteria
                    // Criteria 1: modulus stabilize within 5% at all Gaussian points (less strict criteria only checks the center Gaussian point)
                    VectorXd error = (modulus - modulus_new).array() / modulus_old.array();
                    error = error.array().abs(); // or error.cwiseAbs()
                    if (g == 4 && error(0) > 0.05 && error(1) > 0.05 && error(2) > 0.05) // tutu uses modulus_old, but I want to use modulus
                        convergence = false;
                    // Criteraia 2: Accumulative modulus error within 0.2%
                    error = error.array().square();
                    modulus_old = modulus_old.array().square();
                    if (g == 4/*true*/) { // less strict convergence criteria
                        sumError += error.sum();
                        sumModulus += modulus_old.sum(); // tutu uses modulus_old, but I want to use modulus
                    }

                }
            }

        }

    }
    // std::cout << "Sum Error: " << sumError / sumModulus << std::endl;
    //std::cout << "Modulus Element No.1: " << mesh.elementArray()[1]->modulusAtGaussPt(1) << std::endl;
    return (sumError / sumModulus < 0.002 && convergence) ? true : false;

}

bool Nonlinear::noTensionIteration()
{
    bool convergence = true;

    Element* curr;
    int numNodes; // number of nodes belong to the element
    int numGaussianPt; // number of Gaussian points of the element
    for (int i = 0; i < mesh.elementCount(); i++) {
        curr = mesh.elementArray()[i];
        Material* material = curr->material();
        if (material->noTension) { // compute stress for no tension granular elements only, skip all HMA/subgrade ones
            const VectorXi & nodeList = curr->getNodeList();
            numNodes = curr->getSize();
            numGaussianPt = (int)curr->shape()->gaussianPt().size();

            // Assemble the nodal displacement vector for this element
            VectorXd nodeDisp(2 * numNodes);
            for (int j = 0; j < numNodes; j++) {
                nodeDisp(2 * j) = nodalDisp(2 * nodeList(j));
                nodeDisp(2 * j + 1) = nodalDisp(2 * nodeList(j) + 1);
            }

            // Step 1: Compute stress at gaussian points based on the solved nodal displacement from last iteration
            // Step 2: Compute principal stress and rotation angle based on the stress (sigma r, sigma theta, sigma z, tau rz), and filter out the tensile components based on the limiting tensile strength critera 0.1.
            // Step 3: Counteract the global load vector based on the tensionForce and output the boolean convergence.
            MatrixXd tension(4, numGaussianPt);
            for (int g = 0; g < numGaussianPt; g++) {
                MatrixXd B = curr->BMatrix(curr->shape()->gaussianPt(g));
                VectorXd strain = B * nodeDisp; // e = B * u
                double modulus = (curr->modulusAtGaussPt)(g);
                VectorXd modulus_vec(1);
                modulus_vec << modulus;
                VectorXd stress = material->EMatrix(modulus_vec) * (strain - curr->thermalStrain()); // sigma = E_(i-1) * (e - e0), note that the M and E are both from previous iteration

                // Compute principal stress and rotation angle
                // p1 = (s1 + s3) / 2 + radius;
                // p3 = (s1 + s3) / 2 - radius;
                // p2 = s2;
                // tan(2*theta) = -2 * tau / (s1 - s3)
                double radius = std::sqrt( (stress(0) - stress(2)) * (stress(0) - stress(2)) / 4 + stress(3) * stress(3) ); // sqrt{ [(s1 - s3)/2]^2 + tau^2 }
                double sigma1 = (stress(0) + stress(2)) / 2 + radius;
                double sigma2 = stress(1);
                double sigma3 = (stress(0) + stress(2)) / 2 - radius;
                double theta = std::atan2(-2 * stress(3), stress(0) - stress(2)) / 2;
                // Or use the eigen approach:
                // VectorXd principal = principalStress(stress);
                // double sigma1 = principal(2);
                // double sigma2 = principal(1);
                // double sigma3 = principal(0);
                // double theta = std::atan2(-2 * stress(3), stress(0) - stress(2)) / 2;

                // In our FEM routine, + is tension, - is compression
                double limit = 0.0;
                if (g == 4 && (sigma1 > limit || sigma2 > limit || sigma3 > limit)) // only check at center Gaussian point, otherwise the solution will diverge
                    convergence = false;

                double t1 = sigma1 > limit ? sigma1 : 0;
                double t2 = sigma2 > limit ? sigma2 : 0;
                double t3 = sigma3 > limit ? sigma3 : 0;
                // Principal stress -> Polar stress -> Counteract force
                double sigma_r = (t1 + t3) / 2 + (t1 - t3) * std::cos(2 * theta) / 2;
                double sigma_t = t2;
                double sigma_z = (t1 + t3) / 2 - (t1 - t3) * std::cos(2 * theta) / 2;
                double tau_rz = (t1 - t3) * std::sin(2 * theta) / 2;
                VectorXd tensionStress(4);
                tensionStress << sigma_r, sigma_t, sigma_z, tau_rz;
                tension.col(g) = tensionStress;

                if (i == 107 && g == 1) { // the granular element at centerline
                    for (int mm = 0; mm < numNodes; mm++)
                        std::cout << " NodalForce "<< mm <<": "<< nodalForce(2*nodeList(mm)) <<" "<<nodalForce(2*nodeList(mm)+1) << std::endl;
                }
            }
            VectorXd tensionForce = curr->computeTensionForce(tension);

            // Assemble the element tension force to global nodal force
            for (int k = 0; k < numNodes; k++) {
                nodalForce(2 * nodeList(k)) += tensionForce(2 * k);
                nodalForce(2 * nodeList(k) + 1) += tensionForce(2 * k + 1);
            }
        }
    }
    return convergence;
}

VectorXd Nonlinear::principalStress(const VectorXd & stress) const
{
    // In our coordinates, vertical stress: -:compression +:tension; Horizontal stress: -:compression +:tension
    MatrixXd tensor(3,3);
    tensor << stress(0), 0, stress(3),
              0, stress(1), 0,
              stress(3), 0, stress(2);
    SelfAdjointEigenSolver<MatrixXd> es(tensor, EigenvaluesOnly);
    return es.eigenvalues();

    // Tutu's approach
    // VectorXd result(3);
    // double radius = std::sqrt( (stress(0) - stress(2)) * (stress(0) - stress(2)) / 4 + stress(3) * stress(3) ); // sqrt{ [(s1 - s3)/2]^2 + tau^2 }
    // double sigma1 = (stress(0) + stress(2)) / 2 + radius;
    // double sigma2 = stress(1);
    // double sigma3 = (stress(0) + stress(2)) / 2 - radius;
    // if (sigma2 < sigma3)
    //     std::swap(sigma2, sigma3);
    // result << sigma3, sigma2, sigma1;
    // return result;
}
