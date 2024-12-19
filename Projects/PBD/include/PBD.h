#ifndef PBD_H
#define PBD_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <random>
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

    // data structure for collision constraints
    struct CollisionConstraint
    {
        bool inBoat=false; //is the collision with the boat?

        int fromBelow;

        int qIdx; // vertex

        IV f; // face

        TV n = {0, 0, 0}; // normal
    };

    // data structure for collision constraints with static objects (in our scene, the boat)
    struct CollisionConstraintStatic
    {

        int pIdx; // vertex

        TV q; // face

        TV n = {0, 0, 0}; // normal
    };

public:
    Eigen::MatrixXd boatV;
    Eigen::MatrixXi boatF;

    MatrixXT currentBoatV;        // boat vertex positions, 'x' in paper, (V x 3)

    // vector with 1/mass of each vertex, (V x 1)
    VectorXT boatW;

    // vector with mass of each vertex, (V x 1)
    VectorXT boatM;

    // positions, (V x 3)
    MatrixXT boatP;

    // velocities, (V x 3)
//    MatrixXT v;

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

    // vector with 1/mass of each vertex, (V x 1)
    VectorXT w;

    // vector with mass of each vertex, (V x 1)
    VectorXT m;

    // vector with Lagrange multipliers, (C x 1)
    VectorXT lambdas;

    // map containing constraint name -> constraint idx
    std::unordered_map<std::string, int> constraint_idx = {
        std::make_pair("stretch", 0),
        std::make_pair("bend", 1)
    };

    // time step
    T dt = 0.01;

    // cloth density [kg/m^2]
    T rho = 1.0;

    // cloth thickness [m]
    T h = 0.001;

    TV testPoint;

//    //minimum distance from rigid body
//    T minDist=0

    T maxEdgeLength=0;

    // friction coefficient
    T mu = 0.5;

    // stiffness parameters
    T k_stretch = 0.5;
    T k_bend = 0.25;

    // compliance parameters
    T alpha_stretch = 0.075;
    T alpha_bend = 4.0;

    //damping parameter
    T k_damping = 0.01;

    // gravitational acceleration, (3 x 1)
    TV g = {0.0, 0.0, -9.81};

    // number of iterations to solve the constraints
    size_t numIterations = 30;

    // number of steps (outer loop)
    size_t nSteps=10000;

    // collision constraints
    std::vector<CollisionConstraint> collisionConstraintsList;

    // static collision constraints
    std::vector<CollisionConstraintStatic> collisionConstraintsStaticList;

    // position contraints
    // indices, (V' x 1)
    VectorXi posConstraintsIdxs;
    // vertex positions, (V' x 3)
    MatrixXT posConstraintsV;

    TV constantVelocity={0.5,0,0};

public:
    bool stretchingConstraintsActivated = true;
    bool bendingConstraintsActivated = true;
    bool collisionConstraintsActivated = true;
    bool positionConstraintsActivated = true;
    bool useSpatialHashing = true;
    bool floorCollision = true;
    bool fakeWindActivated = true;
    bool useXPBD = true;

private:
    int seed = 43854397;
    std::mt19937 gen = std::mt19937(seed);
    std::uniform_int_distribution<int> dist10;
    std::uniform_int_distribution<int> distV;
    std::uniform_int_distribution<int> distF;


    int prime1=73856093, prime2=19349663, prime3=83492791;
    float epsilon=0.001f;
    std::vector<std::vector<int>> incidentFaces;
    std::vector<std::vector<int>> adjList;
    int hash(int i, int j, int k, int n);
    int hash(const PBD::TV point, T l, int n, const TV& minCoord);

    static int isAbove(const PBD::TV& q, const PBD::TV& p1, const PBD::TV& p2, const PBD::TV& p3) ;

    static bool pointIntersectsTriangle(const PBD::TV& q, const PBD::TV& p1, const PBD::TV& p2, const PBD::TV& p3, T thickness) ;

    T closestToPointInTriangle(const PBD::TV& q, const PBD::TV& p1, const PBD::TV& p2, const PBD::TV& p3, PBD::TV& closestPoint);

    bool edgesAreClose(const PBD::TV& x1, const PBD::TV& x2, const PBD::TV& x3, const PBD::TV& x4) const;

    void hashVertices(std::vector<std::vector<int>>& hashTable, T boxSize, TV& minCoord);

    void spatialHashing();
    void spatialHashingStatic();
public:
    void initializeFromFile(const std::string& filename);

    void stretchingConstraintsXPBD();
    void stretchingConstraints(int solver_it);

    void bendingConstraintsXPBD();
    void bendingConstraints(int solver_it);

    void generateCollisionConstraints();

    void generateCollisionConstraintsStatic();

    void collisionConstraints();

    void collisionConstraintsStatic();

    void applyFriction();

    void positionConstraints();

    void projectConstraints(int solver_it);

    void dampVelocities(T kDamping);

    void applyFakeWind();

    bool advanceOneStep(int step);

public:
    PBD() {}
    ~PBD() {}

};

#endif
