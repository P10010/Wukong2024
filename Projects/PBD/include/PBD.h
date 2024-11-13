#ifndef PBD_H
#define PBD_H


#include <utility> 
#include <iostream>
#include <fstream>
#include <Eigen/Geometry>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <tbb/tbb.h>
#include <unordered_map>
#include <iomanip>


#include "VecMatDef.h"

class PBD
{
public:
	using VectorXT = Matrix<T, Eigen::Dynamic, 1>;
	using MatrixXT = Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
	using VectorXi = Vector<int, Eigen::Dynamic>;
	using MatrixXi = Matrix<int, Eigen::Dynamic, Eigen::Dynamic>;
	using VtxList = std::vector<int>;
	using StiffnessMatrix = Eigen::SparseMatrix<T>;
	using Entry = Eigen::Triplet<T>;
	using TV = Vector<T, 3>;
	using TV2 = Vector<T, 2>;
	using TM2 = Matrix<T, 2, 2>;
	using TV3 = Vector<T, 3>;
	using IV = Vector<int, 3>;
	using IV2 = Vector<int, 2>;
	using TM = Matrix<T, 3, 3>;

public:
	VectorXi faces;
	VectorXT atRest;
    VectorXT currentV;

    //velocities
    VectorXT v;

	//positions
	VectorXT p, x;

    //vector with 1/mass of each vertex
    VectorXT w;

    //time step
    T dt = 0.01;

    //gravitational acceleration
    T g=9.81;

	//number of iterations to solve the constraints
	size_t numIterations = 10;

public:
	void initializeFromFile(const std::string& filename);

	void projectConstraints();

    bool advanceOneStep(int step);
public:
	PBD() {} 
	~PBD() {} 
};


#endif
