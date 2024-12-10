#ifndef APP_H
#define APP_H

#include "igl/readOBJ.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "Util.h"
#include "polyscope/types.h"
#include "polyscope/view.h"

template <class Simulation>
class App
{
public:
    Simulation& simulation;
    int static_solve_step = 0;

    polyscope::SurfaceMesh* psMesh;

    T t = 0;

    bool animate_modes = false;
    bool run_sim = false;
    int modes = 0;
    Eigen::MatrixXd eigen_vectors;
    Eigen::VectorXd eigen_values;

    Eigen::MatrixXd meshV;
    Eigen::MatrixXi meshF;

    Eigen::MatrixXd boatV;
    Eigen::MatrixXi boatF;

public:
    void initializeScene()
    {

        polyscope::options::autocenterStructures = false;
        polyscope::options::autoscaleStructures = false;
        polyscope::view::upDir = polyscope::UpDir::ZUp;
        polyscope::view::frontDir = polyscope::FrontDir::YFront;
        polyscope::view::windowWidth = 3000;
        polyscope::view::windowHeight = 2000;
        polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
        polyscope::options::groundPlaneHeightFactor = 0.6;
        polyscope::options::shadowDarkness = 0.4;
        // Initialize polyscope
        polyscope::init();
        // TODO: why was this necessary?
        // vectorToIGLMatrix<T, 3>(simulation.currentV, meshV);
        // vectorToIGLMatrix<int, 3>(simulation.faces, meshF);
        meshV = simulation.currentV;
        meshF = simulation.faces;
        // TODO: =========================
        psMesh = polyscope::registerSurfaceMesh("sail", meshV, meshF);
        psMesh->setSmoothShade(false);
        psMesh->setSurfaceColor(glm::vec3(252. / 255., 247. / 255., 216. / 255.));
        psMesh->setEdgeWidth(1.0);

        // TODO: change to PBD
        igl::readOBJ("../../../Projects/PBD/data/sailboat/parts/boat.obj", boatV, boatF);
        polyscope::SurfaceMesh* boatMesh = polyscope::registerSurfaceMesh("boat", boatV, boatF);
        boatMesh->setSmoothShade(false);
        boatMesh->setSurfaceColor(glm::vec3(161. / 255., 77. / 255., 34. / 255.));
        boatMesh->setEdgeWidth(0.0);

        polyscope::state::userCallback = [&]() { sceneCallback(); };
    }
    void sceneCallback()
    {
        if (ImGui::Button("Toggle Simulation"))
        {
          static_solve_step=0;
            run_sim = !run_sim;
        }

        // Reset Simulation
        if (ImGui::Button("Reset Simulation"))
        {
            simulation.initializeFromFile(simulation.scene);
            initializeScene();
            static_solve_step=0;
            run_sim = false;
        }

        // Change stretching, bend and damping parameters
        const double min = 0.0;
        const double max = 1.0;
        ImGui::SliderScalar("Stretching compliance", ImGuiDataType_Double, &simulation.alpha_stretch, &min, &max, "%.3f");
        ImGui::SliderScalar("Bend compliance", ImGuiDataType_Double, &simulation.alpha_bend, &min, &max, "%.3f");
        ImGui::SliderScalar("Damping", ImGuiDataType_Double, &simulation.k_damping, &min, &max, "%.2f");

        // Change rho
        const double min_rho = 0.0;
        const double max_rho = 10.0;
        ImGui::SliderScalar("Density", ImGuiDataType_Double, &simulation.rho, &min_rho, &max_rho, "%.2f");

        // Change thickness h
        const double min_h = 0.0001;
        const double max_h = 0.01;
        ImGui::SliderScalar("Thickness", ImGuiDataType_Double, &simulation.h, &min_h, &max_h, "%.5f");

        // Change number of iterations
        const size_t min_iter = 0;
        const size_t max_iter = 100;
        ImGui::SliderScalar("Solver iterations", ImGuiDataType_U64, &simulation.numIterations, &min_iter, &max_iter, "%d");

        // Change time step
        const double dt_step = 0.001, dt_stepfast = 0.01;
        ImGui::InputDouble("Time step [s]", &simulation.dt, dt_step, dt_stepfast);

        // Change number of steps
        const size_t min_steps = 1;
        const size_t max_steps = 10000;
        ImGui::SliderScalar("Steps", ImGuiDataType_U64, &simulation.nSteps, &min_steps, &max_steps, "%d");

        // Activate/Deactivate constraints
        ImGui::Checkbox("Position Constraints", &simulation.positionConstraintsActivated);
        ImGui::Checkbox("Stretching Constraints", &simulation.stretchingConstraintsActivated);
        ImGui::Checkbox("Bending Constraints", &simulation.bendingConstraintsActivated);
        ImGui::Checkbox("Collision Constraints", &simulation.collisionConstraintsActivated);
        ImGui::Checkbox("Spatial Hashing", &simulation.useSpatialHashing);
        ImGui::Checkbox("Floor Collision", &simulation.floorCollision);

        if (!animate_modes && run_sim)
        {
            bool finished = simulation.advanceOneStep(static_solve_step++);
            Eigen::MatrixXd meshV = simulation.currentV;
            // TODO: why was this necessary?
            // vectorToIGLMatrix<T, 3>(simulation.currentV, meshV);
            psMesh->updateVertexPositions(meshV);
            if (finished)
                run_sim = false;
        }
    }
    void run() { polyscope::show(); }

public:
    App(Simulation& sim) : simulation(sim) {}
    ~App() {}
};

#endif
