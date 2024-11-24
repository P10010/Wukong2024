#ifndef APP_H
#define APP_H

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

public:
    void initializeScene()
    {
        polyscope::options::autocenterStructures = true;
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
        psMesh = polyscope::registerSurfaceMesh("surface mesh", meshV, meshF);
        psMesh->setSmoothShade(false);
        psMesh->setSurfaceColor(glm::vec3(0.255, 0.514, 0.996));
        psMesh->setEdgeWidth(1.0);
        polyscope::state::userCallback = [&]() { sceneCallback(); };
    }
    void sceneCallback()
    {
        if (ImGui::Button("Toggle Simulation"))
        {
            run_sim = !run_sim;
        }

        // Reset Simulation
        if (ImGui::Button("Reset Simulation"))
        {
            simulation.initializeFromFile(simulation.scene);
            initializeScene();
            run_sim = false;
        }

        // Change stretching and bend parameters
        const double min = 0.0;
        const double max = 1.0;
        ImGui::SliderScalar("Stretching", ImGuiDataType_Double, &simulation.k_stretch, &min, &max, "%.2f");
        ImGui::SliderScalar("Bend", ImGuiDataType_Double, &simulation.k_bend, &min, &max, "%.2f");

        // Change rho
        const double min_rho = 0.0;
        const double max_rho = 10.0;
        ImGui::SliderScalar("Density", ImGuiDataType_Double, &simulation.rho, &min_rho, &max_rho, "%.2f");

        // Change number of iterations
        const size_t min_iter = 0;
        const size_t max_iter = 100;
        ImGui::SliderScalar("Collision iterations", ImGuiDataType_U64, &simulation.numIterations, &min_iter, &max_iter, "%d");

        // Change time step
        const double dt_step=0.001,dt_stepfast=0.01;
        ImGui::InputDouble("Time step [s]", &simulation.dt, dt_step, dt_stepfast);

        // Change number of steps
        const size_t min_steps = 50;
        const size_t max_steps = 10000;
        ImGui::SliderScalar("Steps", ImGuiDataType_U64, &simulation.nSteps, &min_steps, &max_steps, "%d");

        // Activate/Deactivate constraints
        ImGui::Checkbox("Stretching Constraints", &simulation.stretchingConstraintsActivated);
        ImGui::Checkbox("Bending Constraints", &simulation.bendingConstraintsActivated);
        ImGui::Checkbox("Collision Constraints", &simulation.collisionConstraintsActivated);
        ImGui::Checkbox("Position Constraints", &simulation.positionConstraintsActivated);

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
