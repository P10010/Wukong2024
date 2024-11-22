#ifndef PBD_H
#define PBD_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <tbb/tbb.h>
#include <unordered_map>
#include <utility>

#include "VecMatDef.h"

class PBD
{
public:
    using VectorXT = VMD::Matrix<T, Eigen::Dynamic, 1>;
    using MatrixXT = VMD::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using VectorXi = VMD::Vector<int, Eigen::Dynamic>;
    using MatrixXi = VMD::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>;
    using VtxList = std::vector<int>;
    using StiffnessMatrix = Eigen::SparseMatrix<T>;
    using Entry = Eigen::Triplet<T>;
    using TV = VMD::Vector<T, 3>;
    using TV2 = VMD::Vector<T, 2>;
    using TM2 = VMD::Matrix<T, 2, 2>;
    using TV3 = VMD::Vector<T, 3>;
    using IV = VMD::Vector<int, 3>;
    using IV2 = VMD::Vector<int, 2>;
    using TM = VMD::Matrix<T, 3, 3>;

public:
    std::string scene;        // path to the obj file
    MatrixXi faces;           // (F x 3)
    MatrixXi TT;              // triangle-triangle adjacencies, (F x 3)
    MatrixXi TTi;             // triangle-triangle adjacencies, (F x 3)
    MatrixXT dihedral_angles; // dihedral angles between between faces, (F x 3)
    MatrixXi edges;           // (E x 2)
    VectorXT edge_lengths;    // (E x 2)
    MatrixXT atRest;          // initial positions (V x 3)
    MatrixXT currentV;        // vertex positions, 'x' in paper, (V x 3)

    // velocities, (V x 3)
    MatrixXT v;

    // positions, (V x 3)
    MatrixXT p;

    // position contraints
    // indices, (V' x 1)
    VectorXi posConstraintsIdxs;
    // vertex positions, (V' x 3)
    MatrixXT posConstraintsV;

    // vector with 1/mass of each vertex, (V x 1)
    VectorXT w;

    // time step
    T dt = 0.01;

    // cloth density [kg/m^2]
    T rho = 1.0;

    // stiffness parameters
    T k_stretch = 0.5;
    T k_bend = 0.25;

    // gravitational acceleration, (3 x 1)
    TV g = {0.0, 0.0, -9.81};

    // number of iterations to solve the constraints
    size_t numIterations = 10;


public:
    bool stretchingConstraintsActivated = true;
    bool bendingConstraintsActivated = true;
    bool positionConstraintsActivated = true;

public:
    void initializeFromFile(const std::string& filename);

    void stretchingConstraints(int solver_it);

    void bendingConstraints(int solver_it);

    void positionConstraints();

    void projectConstraints(int solver_it);

    bool advanceOneStep(int step);

public:
    PBD() {}
    ~PBD() {}
};

#endif
