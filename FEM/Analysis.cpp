/**
 * @file Analysis.cpp
 * Implementation of Analysis class.
 *
 * @author Haohang Huang
 * @date Feburary 26, 2018
 */

#define _USE_MATH_DEFINES
#include <cmath>
#include "Analysis.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <unordered_map>

// Analysis::Analysis()
// {
// }

Analysis::Analysis(Mesh & meshInfo)
  :
  mesh(meshInfo), globalStiffness(2 * mesh.nodeCount(), 2 * mesh.nodeCount()),
  nodalDisp(VectorXd::Zero(2 * mesh.nodeCount())), nodalForce(VectorXd::Zero(2 * mesh.nodeCount())),
  nodalStrain(MatrixXd::Zero(mesh.nodeCount(), 4)), nodalStress(MatrixXd::Zero(mesh.nodeCount(), 4))
{
}

Analysis::~Analysis()
{
}

void Analysis::assembleStiffness()
{
    // Reserve the space by estimating the maximum number of non-zero elements
    // in the sparse matrix. This should actually be No. of elements * 2n * 2n
    // where n is the size of the element type with the most nodes. Here n = 8.
    // unsigned long int maxNum = mesh.elementCount() * 16 * 16;
    // globalStiffness.reserve(maxNum);
    //
    // // Assemble global matrix from local matrix of each element
    // Element* curr;
    // for (int i = 0; i < mesh.elementCount(); i++) {
    //     curr = mesh.elementArray()[i];
    //     int size = curr->getSize();// element type, for Q4 element, size=4; for Q8, size=8, etc
    //     const VectorXi & nodeList = curr->getNodeList();// the index of nodes belong to this element, e.g., for element8, it will give you a vector contain (10,11,15,14), use this to locate row & column in globalStiffness matrix
    //     const MatrixXd & localStiffness = curr->localStiffness();
    //
    //     // Traverse each node and assemble the values to global stiffness matrix and global force vector
    //     for (int j = 0; j < size; j++){
    //
    //         for (int k = 0; k < size; k++){
    //             globalStiffness.coeffRef(2 * nodeList(j), 2 * nodeList(k)) += localStiffness(2 * j , 2 * k);
    //             globalStiffness.coeffRef(2 * nodeList(j), 2 * nodeList(k) + 1) += localStiffness(2 * j , 2 * k + 1);
    //             globalStiffness.coeffRef(2 * nodeList(j) + 1, 2 * nodeList(k)) += localStiffness(2 * j + 1, 2 * k);
    //             globalStiffness.coeffRef(2 * nodeList(j) + 1, 2 * nodeList(k) + 1) += localStiffness(2 * j + 1, 2 * k + 1);
    //         }
    //
    //         // Also assemble the body force and temperature load vector
    //         const VectorXd & forceVec = curr->nodalForce();
    //         nodalForce(2 * nodeList(j)) += forceVec(2 * j);
    //         nodalForce(2 * nodeList(j) + 1) += forceVec(2 * j + 1);
    //
    //     }
    //
    // }
    // globalStiffness.makeCompressed();

    // Sparse matrix operation notes:
    //     m.setZero() to remove all non-zero coefficients
    //     m.rows() to get number of rows
    //     m.cols() to get number of columns
    //     m.coeffRef(i,j) = k to set value to the element already exists;
    //     m.insert(i,j) = k to set value to the element does not already exist;
    //
    // Efficiently create sparse matrix from triplet list/vector (i,j,value):
    //     typedef Eigen::Triplet<double> T;
    //     std::vector<T> tripletList;
    //     int estimation_of_entries;
    //     tripletList.reserve(estimation_of_entries);
    //     tripletList.push_back(T(i,j,value)); // or emplace_back()
    //     globalStiffness.setFromTriplets(tripletList.begin(), tripletList.end());
    // *Limitation: once initialized, the value cannot be modified

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    unsigned long int maxNum = mesh.elementCount() * 16 * 16;
    tripletList.reserve(maxNum);

    const std::vector<int> & DOFList = mesh.boundaryNodeList;
    const std::vector<double> & boundaryValue = mesh.boundaryValue;
    std::unordered_map<unsigned long int, double> boundaryHash;
    boundaryHash.reserve(DOFList.size());
    for (unsigned i = 0; i < DOFList.size(); i++)
        boundaryHash[DOFList[i]] = boundaryValue[i];

    // Assemble global matrix from local matrix of each element, meanwhile modify
    // the stiffness matrix based on boundary condition and adjust the force vector
    // accordingly
    Element* curr;
    for (int i = 0; i < mesh.elementCount(); i++) {
        curr = mesh.elementArray()[i];
        int size = curr->getSize();// element type, for Q4 element, size=4; for Q8, size=8, etc
        const VectorXi & nodeList = curr->getNodeList();// the index of nodes belong to this element, e.g., for element8, it will give you a vector contain (10,11,15,14), use this to locate row & column in globalStiffness matrix

        // Bootstrap the computation of local stiffness matrix and force vector. After calling this function, the member variables are all computed
        // @BUG (solved) previous this bootstrap step is in the ctor of derived class ElementQ8, so in the nonlinear analysis, the localStiffness and body & temp force are only computed once at the beginning!
        curr->computeStiffnessAndForce();
        // if (i == 0)
        //     std::cout << "local stiff: " << curr->localStiffness() << std::endl;
        const MatrixXd & localStiffness = curr->localStiffness();

        // Traverse each node and assemble the values to global stiffness matrix and global force vector
        for (int j = 0; j < size; j++){
            for (int k = 0; k < size; k++) {
                // Since the j,k are the node index, actually for every (j,k) we
                // are accessing a 2x2 block in the local stiffness matrix:
                // (2j,2k), (2j+1,2k), (2j+1,2k+1), (2j,2k+1)
                // For applying the boundary condition, we subtract the column
                // multiplied by the boundar value, and then cross out the column
                // and row at the fixed DOF. So rows are useless and we don't
                // assign into sparse matrix. Columns are subtracted incrementally
                // from the force vector. A hash table is used to check the DOF
                // location in constant time.
                // Think from a column-wise perspective
                // First column
                if (boundaryHash.find(2 * nodeList(k)) != boundaryHash.end()) { // if on the crossed-out column, subtract it from force vector
                    if (j != k) // if at the crossing (j = k), assign as 1 at the end. Only (2j,2k) can possibly be on the diagonal
                        nodalForce(2 * nodeList(j)) -= localStiffness(2 * j , 2 * k) * boundaryHash[2 * nodeList(k)];
                    nodalForce(2 * nodeList(j) + 1) -= localStiffness(2 * j + 1, 2 * k) * boundaryHash[2 * nodeList(k)];
                }
                else {
                    if (boundaryHash.find(2 * nodeList(j)) == boundaryHash.end()) // if on the crossed-out row, do nothing; otherwise add it to the sparse K
                        tripletList.push_back(T(2 * nodeList(j), 2 * nodeList(k), localStiffness(2 * j , 2 * k)));
                    if (boundaryHash.find(2 * nodeList(j) + 1) == boundaryHash.end())
                        tripletList.push_back(T(2 * nodeList(j) + 1, 2 * nodeList(k), localStiffness(2 * j + 1 , 2 * k)));
                }

                // Second column
                if (boundaryHash.find(2 * nodeList(k) + 1) != boundaryHash.end()) {
                    if (j != k) // if at the crossing (j = k), assign as 1 at the end. Only (2j+1,2k+1) can possibly be on the diagonal
                        nodalForce(2 * nodeList(j) + 1) -= localStiffness(2 * j + 1, 2 * k + 1) * boundaryHash[2 * nodeList(k) + 1];
                    nodalForce(2 * nodeList(j)) -= localStiffness(2 * j , 2 * k + 1) * boundaryHash[2 * nodeList(k) + 1];
                }
                else {
                    if (boundaryHash.find(2 * nodeList(j)) == boundaryHash.end()) // if on the crossed-out row, do nothing; otherwise add it to the sparse K
                        tripletList.push_back(T(2 * nodeList(j), 2 * nodeList(k) + 1, localStiffness(2 * j , 2 * k + 1)));
                    if (boundaryHash.find(2 * nodeList(j) + 1) == boundaryHash.end())
                        tripletList.push_back(T(2 * nodeList(j) + 1, 2 * nodeList(k) + 1, localStiffness(2 * j + 1 , 2 * k + 1)));
                }

            }

            // Also assemble the body force and temperature load vector element-wise
            // meanwhile finalize the force value at boundary locations (the point
            // and edge load are applied in the steps above, so here after we apply the body
            // force and temperature load, we should finalize the boundary locations
            // in the global force vector to be boundary value to get 1 x U = U effect
            const VectorXd & forceVec = curr->nodalForce();
            if (boundaryHash.find(2 * nodeList(j)) == boundaryHash.end())
                nodalForce(2 * nodeList(j)) += forceVec(2 * j);
            else
                nodalForce(2 * nodeList(j)) = boundaryHash[2 * nodeList(j)];
            if (boundaryHash.find(2 * nodeList(j) + 1) == boundaryHash.end())
                nodalForce(2 * nodeList(j) + 1) += forceVec(2 * j + 1);
            else
                nodalForce(2 * nodeList(j) + 1) = boundaryHash[2 * nodeList(j) + 1];

        }

    }

    // Assign the crossing of boundary locations to 1, so the corresponding
    // location in U and F are both assigned the boundary value, we have 1 x U = U
    // and can always have the exact answer at the boundary node
    for (unsigned i = 0; i < DOFList.size(); i++)
        tripletList.push_back(T(DOFList[i], DOFList[i], 1));

    // Write into sparse matrix
    globalStiffness.setFromTriplets(tripletList.begin(), tripletList.end());

    // Assign the crossing of boundary locations to 1
    // You can also skip every crossing in the following steps, and here set a 1 to the triplet
    // for (unsigned i = 0; i < DOFList.size(); i++)
    //     globalStiffness.coeffRef(DOFList[i], DOFList[i]) = 1;

    globalStiffness.makeCompressed();

    // Free the memory
    std::vector<T>().swap(tripletList);
}

void Analysis::applyForce()
{
    // @BUG (solved) previous miss this initialization step, therefore in the nonlinear analysis the force will accumulate in every iteration and blow up!!!
    // Other variable such as nodelDisp, globalStiffness, nodalStrain, nodalStress will be rewritten every time, so doesn't matter
    nodalForce = VectorXd::Zero(2 * mesh.nodeCount());

    // Apply point load
    for (unsigned i = 0; i < mesh.loadNodeList.size(); i++)
        nodalForce(mesh.loadNodeList[i]) += 2 * M_PI * mesh.loadValue[i]; // *2 M_PI due to the axisymmetric property

    // Apply edge load
    Element* curr;
    // Traverse all elements with edge load
    for (unsigned i = 0; i < mesh.loadElementList.size(); i++) {
        curr = mesh.getElement(mesh.loadElementList[i]);
        const VectorXi & nodeList = curr->getNodeList();
        int elementType = curr->getSize();
        // The shape function for each edge in an isoparametric element is the same
        const std::vector<double> & gaussianPoint = curr->shape()->edgeGaussianPt(); // length 3 vector
        const std::vector<double> & gaussianWeight = curr->shape()->edgeGaussianWt(); // length 3 vector
        int numGaussianPts = (int)gaussianPoint.size();

        // The edge load information
        const std::vector<int> & edges = mesh.loadEdgeList[i];
        const std::vector<double> & loads = mesh.edgeLoadValue[i];

        VectorXd forceVec(2 * elementType);
        forceVec.setZero();
        // Traverse all loaded edges in this element
        for (unsigned j = 0; j < edges.size(); j++) {

            Vector2d load(loads[2 * j], loads[2 * j + 1]); // 2x1 load vector

            // The global coordinates of the edge nodes
            const std::vector<int> & edgeNodes = curr->shape()->edge(edges[j]);
            int numEdgeNodes = (int)edgeNodes.size();
            MatrixXd nodeCoord (2, numEdgeNodes);
            for (int n = 0; n < numEdgeNodes; n++)
                nodeCoord.col(n) = mesh.getNode(nodeList(edgeNodes[n]))->getGlobalCoord();
                // (Solved) @BUG here, previous wrote mesh.getNode(edgeNodes[n]), should map to the node index in the mesh

            // Integration at gaussian points
            // sum: 2PI * N^T * F * |J| * r * W(i)
            VectorXd result(2 * numEdgeNodes); result.setZero(); // 6x1 vector
            for (int g = 0; g < numGaussianPts; g++) {
                MatrixXd N = curr->shape()->edgeFunctionMat(g); // 2x6 shape matrix
                VectorXd globalDeriv = nodeCoord * curr->shape()->edgeFunctionDeriv(g); // 2x3 matrix * 3x1 vector = 2x1 vector [dr; dz]
                double dr = globalDeriv(0);
                double dz = globalDeriv(1);
                double jacobianDet = std::sqrt(dr * dr + dz * dz);
                double radius = nodeCoord.row(0) * curr->shape()->edgeFunctionVec(g);
                result += 2 * M_PI * N.transpose() * load * jacobianDet * radius * gaussianWeight[g];
            }

            // Assemble edge force vector to the element force vector
            for (int n = 0; n < numEdgeNodes; n++) {
                forceVec(2 * edgeNodes[n]) += result(2 * n);
                forceVec(2 * edgeNodes[n] + 1) += result(2 * n + 1);
            }

        } // After traverse all loaded edges of this element

        // Assemble element force vector to the global force vector
        for (int k = 0; k < elementType; k++) {
            nodalForce(2 * nodeList(k)) += forceVec(2 * k);
            nodalForce(2 * nodeList(k) + 1) += forceVec(2 * k + 1);
        }

    }

}


// Old version of applying boundary conditon (very slow because we are modifying the whole sparse matrix)
// void Analysis::boundaryCondition()
// {
//     const std::vector<int> & DOFList = mesh.boundaryNodeList;
//     const std::vector<double> & boundaryValue = mesh.boundaryValue;
//     for (unsigned i = 0; i < DOFList.size(); i++) {
//         // Modify stiffness matrix
//         double temp = globalStiffness.coeffRef(DOFList[i], DOFList[i]);
//         for (int j = 0; j < globalStiffness.cols(); j++) {
//             globalStiffness.coeffRef(DOFList[i], j) = 0;
//             globalStiffness.coeffRef(j, DOFList[i]) = 0;
//         } // the sparse matrix .row() or .col() is read-only, we can only do element-wise assignment
//         globalStiffness.coeffRef(DOFList[i], DOFList[i]) = temp;
//
//         // Modify force vector
//         nodalForce -= globalStiffness.col(DOFList[i]) * boundaryValue[i]; // if consider temperature effect : *6.5*-0.000001*30
//         nodalForce(DOFList[i]) = temp * boundaryValue[i];
//     }
// }

void Analysis::computeStrainAndStress()
{
    // Idea:
    // 1. For each element, compute the strain at Gaussian points using nodal
    // displacments by e = Bu.
    // 2. Let the nodal strain be unknowns, from equation we can solve for
    // e(Gaussian) = N e(nodal) --> e(nodal) = N^-1 e(Gaussian)
    // where N^-1 is the pesudo inverse of the N matrix
    // 3. Compute stress from strain by sigma = E e
    // 4. Cumulate strain and stress value at each node, and average at the end
    // In detail of step 2:
    // 9 Gaussian points & 8 nodes, we have e(Gaussian) as a 9-by-4 matrix.
    // N is the stacked shape function [N1 N2 ... N8] at 9 Gaussian points, i.e.
    // a 9-by-8 matrix. And e(nodal) is the stacked strain vector of 8 nodes as
    // a 8-by-4 matrix. Since N is not square, pesudo inverse is used to solve
    // the least square system.

    // Traverse each element
    Element* curr;
    int numNodes; // number of nodes belong to the element
    int numGaussianPt; // number of Gaussian points of the element
    for (int i = 0; i < mesh.elementCount(); i++) {
        curr = mesh.elementArray()[i];
        const VectorXi & nodeList = curr->getNodeList();
        numNodes = curr->getSize();
        numGaussianPt = (int)curr->shape()->gaussianPt().size();


        // Assemble the nodal displacement vector for an element (directly from the solved displacement vector)
        VectorXd nodeDisp(2 * numNodes);
        for (int j = 0; j < numNodes; j++) {
            nodeDisp(2 * j) = nodalDisp(2 * nodeList(j));
            nodeDisp(2 * j + 1) = nodalDisp(2 * nodeList(j) + 1);
        }

        MatrixXd strainAtGaussPt(numGaussianPt, 4); // 4 for axisymmetric problem
        MatrixXd stressAtGaussPt(numGaussianPt, 4);
        MatrixXd shapeAtGaussPt(numGaussianPt, numNodes); // 9x8 matrix for element Q8

        // Compute strain and stress at gaussian points from e = Bu, sigma = Ee
        for (int g = 0; g < numGaussianPt; g++) {
            MatrixXd B = curr->BMatrix(curr->shape()->gaussianPt(g));
            VectorXd e = B * nodeDisp; // e = B * u
            strainAtGaussPt.row(g) = e.transpose();
            VectorXd modulus = (curr->modulusAtGaussPt).row(g); // for nonlinear, this is the stabilized modulus at the Gaussian point; for linear elastic, it's just the constant modulus M
            stressAtGaussPt.row(g) = (curr->EMatrix(modulus) * (e - curr->thermalStrain())).transpose(); // subtract thermal strain, stress = E * (strain - thermal strain)
            shapeAtGaussPt.row(g) = curr->shape()->functionVec(g).transpose();
        }

        // Solve/extrapolate for nodal strain value via a least square linear system using pesudo inverse
        // Previous attempt: solving the system by pesudo inverse, this might have numerical error
        // MatrixXd pesudo = shapeAtGaussPt.completeOrthogonalDecomposition().pseudoInverse(); // pesudo inverse in "Eigen/QR"
        // MatrixXd strainAtNodes = pesudo * strainAtGaussPt; // 8x4 matrix
        // Current solution: use SVD decomposition
        MatrixXd strainAtNodes = shapeAtGaussPt.bdcSvd(ComputeThinU | ComputeThinV).solve(strainAtGaussPt);
        MatrixXd stressAtNodes = shapeAtGaussPt.bdcSvd(ComputeThinU | ComputeThinV).solve(stressAtGaussPt);

        // Notes on LLS system:
        // Several options for solving a linear least squares system:
        // 1. SVD decomposition
        // VectorXd x = A.bdcSvd(ComputeThinU | ComputeThinV).solve(b);
        // 2. QR decomposition
        // VectorXd x = A.colPivHouseholderQr().solve(b);
        // from fast to slow, unstable to stable: householderQr()-->colPivHouseholderQr()-->fullPivHouseholderQr()
        // 3. Normal equations (A^T*A)*x = A^T*b
        // VectorXd x = (A.transpose() * A).ldlt().solve(A.transpose() * b)

        // Set calculated strain and stress value to every node (to be accumulated at each node and averaged later)
        for (int n = 0; n < numNodes; n++) {
            VectorXd strain = strainAtNodes.row(n);
            VectorXd stress = stressAtNodes.row(n);
            mesh.nodeArray()[nodeList(n)]->setStrainAndStress(strain, stress);
        }

    }

}

void Analysis::averageStrainAndStress()
{
    // Average nodal strain and stress, and write the displacement information at the same time!
    for (int i = 0; i < mesh.nodeCount(); i++) {
        mesh.nodeArray()[i]->setDisp(nodalDisp(2 * i), nodalDisp(2 * i + 1));
        nodalStrain.row(i) = mesh.nodeArray()[i]->averageStrain().transpose();
        nodalStress.row(i) = mesh.nodeArray()[i]->averageStress().transpose();
    }
}

void Analysis::printDisp() const
{
    std::cout << "Nodal Displacement: ";
    std::cout << std::endl;
    for (int i = 0; i < mesh.nodeCount(); i++) {
      std::cout << "Node " << mesh.nodeArray()[i]->getIndex() << " : " << mesh.nodeArray()[i]->getDisp().transpose() << std::endl;
    }
    std::cout << std::endl;
}

void Analysis::printStrain() const
{
    std::cout << "Averaged nodal strain: ";
    std::cout << std::endl;
    for (int i = 0; i < mesh.nodeCount(); i++) {
        std::cout << "Node " << i << " : " << nodalStrain.row(i) << std::endl;
    }
    std::cout << std::endl;
    //std::cout << nodalStrain << std::endl; // the whole matrix
}

void Analysis::printStress() const
{
    std::cout << "Averaged nodal stress: ";
    std::cout << std::endl;
    for (int i = 0; i < mesh.nodeCount(); i++) {
        // std::cout << "Node " << i << " : " << nodalStress.row(i) << std::endl;
        VectorXd stress = nodalStress.row(i);
        MatrixXd tensor(3,3);
        tensor << stress(0), 0, stress(3),
                  0, stress(1), 0,
                  stress(3), 0, stress(2);
        SelfAdjointEigenSolver<MatrixXd> es(tensor, EigenvaluesOnly);
        std::cout << "Node principal " << i << " : "<< es.eigenvalues().transpose() << std::endl;
    }
    std::cout << std::endl;
    //std::cout << nodalStress << std::endl; // the whole matrix
}

void Analysis::writeToFile(std::string const & fileName) const
{
    std::ofstream file(fileName);
    std::string END1, END2;

    int temp = std::log10(mesh.nodeCount() - 1);
    if (temp <= 0){
        END1 = "                  ";
        END2 = "        ";
    }
    else if (temp == 1){
        END1 = "                  ";
        END2 = "       ";
    }
    else if (temp == 2){
        END1 = "                  ";
        END2 = "      ";
    }
    else if (temp == 3){
        END1 = "                 ";
        END2 = "      ";
    }
    else if (temp == 4){
        END1 = "                 ";
        END2 = "     ";
    }
    else if (temp == 5){
        END1 = "                 ";
        END2 = "    ";
    }
    else if (temp == 6){
        END1 = "                ";
        END2 = "    ";
    }
    else if (temp == 7){
        END1 = "                ";
        END2 = "   ";
    }
    else if (temp == 8){
        END1 = "                ";
        END2 = " ";
    }
    else {
        END1 = "                ";
        END2 = "";
    }

    file << "\n";
    file << "\n";
    file << "%-------------------------------------------------------------------------------------------------------%\n";
    file << "                  "<< "2D Axisymmetric Finite Element Analysis (Linear Elastic version)\n";
    file << "%-------------------------------------------------------------------------------------------------------%\n";
    file << "\n";
    file << "\n";
    file << "%---------------------------------------------------------------%\n";
    file << "                 "<<"FEM Program Standard Output\n";
    file << "%---------------------------------------------------------------%\n";
    file << "\n";
    file << "                     "<<"Nodal Displacements (in.): \n";
    file << "\n";
    file << "            "<<"  Nodal Index  |      Displacement (R, Z)\n";
    file << "            "<<"--------------------------------------------------\n";
    for (int i = 0; i < mesh.nodeCount(); i++){
        file << END1 <<std::setw(temp + 1) << std::setfill('0')
        << i << END2 <<"|      "<<std::fixed<<std::scientific<<std::setprecision(3)
        <<"("<<mesh.nodeArray()[i]->getDisp().transpose()[0]<<", "<<mesh.nodeArray()[i]->getDisp().transpose()[1]
        <<")\n";
    }
    file << "            "<<"--------------------------------------------------\n";
    file << "\n";
    file << "\n";


    file << "%-------------------------------------------------------------------------------------------------------%\n";
    file << "                                         "<<"FEM Program Standard Output\n";
    file << "%-------------------------------------------------------------------------------------------------------%\n";
    file << "\n";
    file << "                                            "<<"Averaged Nodal Strain (no unit): \n";
    file << "\n";
    file << "            "<<"  Nodal Index  |                     Strain (R, TH, Z, RZ) \n";
    file << "            "<<"--------------------------------------------------------------------------------------\n";
    for (int i = 0; i < mesh.nodeCount(); i++){
        file << END1 <<std::setw(temp + 1) << std::setfill('0')
        << i << END2 <<"|      "<<std::fixed<<std::scientific<<std::setprecision(6)
        <<"("<<nodalStrain.row(i)[0]<<", "<<nodalStrain.row(i)[1]<<", "<<nodalStrain.row(i)[2]<<", "<<nodalStrain.row(i)[3]
        <<")\n";
    }
    file << "            "<<"--------------------------------------------------------------------------------------\n";
    file << "\n";
    file << "\n";

    file << "%-------------------------------------------------------------------------------------------------------%\n";
    file << "                                         "<<"FEM Program Standard Output\n";
    file << "%-------------------------------------------------------------------------------------------------------%\n";
    file << "\n";
    file << "                                            "<<"Averaged Nodal Stress (psi.): \n";
    file << "\n";
    file << "            "<<"  Nodal Index  |                     Stress (R, TH, Z, RZ) \n";
    file << "            "<<"--------------------------------------------------------------------------------------\n";
    for (int i = 0; i < mesh.nodeCount(); i++){
        file << END1 <<std::setw(temp + 1) << std::setfill('0')
        << i << END2 <<"|      "<<std::fixed<<std::scientific<<std::setprecision(6)
        <<"("<<nodalStress.row(i)[0]<<", "<<nodalStress.row(i)[1]<<", "<<nodalStress.row(i)[2]<<", "<<nodalStress.row(i)[3]
        <<")\n";
    }
    file << "            "<<"--------------------------------------------------------------------------------------\n";

    file.close();
}

void Analysis::writeToVTK(std::string const & fileName) const
{
    std::ofstream file(fileName);

    file << "# vtk DataFile Version 2.0\n";
    file << "./vtk\n";
    file << "ASCII\n";
    file << "\n";
    file << "DATASET UNSTRUCTURED_GRID\n";
    file << "POINTS " << mesh.nodeCount() << " double\n";
    for (int i = 0; i < mesh.nodeCount(); i++){
        file << (double)(mesh.nodeArray()[i]->getGlobalCoord())(0) << " " << (double)(mesh.nodeArray()[i]->getGlobalCoord())(1) << " " << "0.0" << "\n";
    }

    file << "\n";

    file << "CELLS " << mesh.elementCount() << " " << (8 + 1) * mesh.elementCount() << "\n";
    for (int i = 0; i < mesh.elementCount(); i++){
        file << 8 << " " << mesh.elementArray()[i]->getNodeList().transpose() << "\n";
    }

    file << "\n";

    file << "CELL_TYPES " << mesh.elementCount() << "\n";
    for (int i = 0; i < mesh.elementCount(); i++){
        file << 23 << "\n"; // quadrilateral
    }

    file << "\n";

    file << "POINT_DATA " << mesh.nodeCount() << "\n";
    // file << "VECTORS " << " Displacement" << " double" << "\n";
    // for (int i = 0; i < nodalDisp.size() / 2; i++){
    //     file << (double)nodalDisp(2 * i) << " " << (double)nodalDisp(2 * i + 1) << " " << 0 << "\n";
    // }

    file << "SCALARS " << "Displacement_R " << "double " << "1" << "\n";
    file << "LOOKUP_TABLE " << "default" << "\n";
    for (int i = 0; i < mesh.nodeCount(); i++){
        file << nodalDisp(2 * i) << "\n";
    }

    file << "SCALARS " << "Displacement_Z " << "double " << "1" << "\n";
    file << "LOOKUP_TABLE " << "default" << "\n";
    for (int i = 0; i < mesh.nodeCount(); i++){
        file << - nodalDisp(2 * i + 1) << "\n";
    }

    file << "SCALARS " << "Stress_R " << "double " << "1" << "\n";
    file << "LOOKUP_TABLE " << "default" << "\n";
    for (int i = 0; i < mesh.nodeCount(); i++)
        file << - nodalStress(i, 0) << "\n";

    file << "SCALARS " << "Stress_Theta " << "double " << "1" << "\n";
    file << "LOOKUP_TABLE " << "default" << "\n";
    for (int i = 0; i < mesh.nodeCount(); i++)
        file << - nodalStress(i, 1) << "\n";

    file << "SCALARS " << "Stress_Z " << "double " << "1" << "\n";
    file << "LOOKUP_TABLE " << "default" << "\n";
    for (int i = 0; i < mesh.nodeCount(); i++)
        file << - nodalStress(i, 2) << "\n";

    file << "SCALARS " << "Stress_Shear " << "double " << "1" << "\n";
    file << "LOOKUP_TABLE " << "default" << "\n";
    for (int i = 0; i < mesh.nodeCount(); i++)
        file << - nodalStress(i, 3) << "\n";

    file << "SCALARS " << "Radial_Distance " << "double " << "1" << "\n";
    file << "LOOKUP_TABLE " << "default" << "\n";
    for (int i = 0; i < mesh.nodeCount(); i++)
        file << (mesh.nodeArray()[i]->getGlobalCoord())(0)<< "\n";

    file << "SCALARS " << "Depth " << "double " << "1" << "\n";
    file << "LOOKUP_TABLE " << "default" << "\n";
    for (int i = 0; i < mesh.nodeCount(); i++)
        file << (mesh.nodeArray()[i]->getGlobalCoord())(1)<< "\n";

    file.close();
}
