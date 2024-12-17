#ifndef APP_H
#define APP_H

#include "Util.h"
#include "igl/readOBJ.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

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

        polyscope::options::autocenterStructures = false;
        polyscope::options::autoscaleStructures = false;
        polyscope::view::upDir = polyscope::UpDir::ZUp;
        polyscope::view::frontDir = polyscope::FrontDir::YFront;
        polyscope::view::windowWidth = 3000;
        polyscope::view::windowHeight = 2000;
        polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
        polyscope::options::groundPlaneHeightFactor = 0.01;
        polyscope::options::shadowDarkness = 0.4;
        polyscope::view::bgColor = {114. / 255., 198. / 255., 243. / 255., 0.};
        polyscope::options::transparencyMode = polyscope::TransparencyMode::Pretty;
        // Initialize polyscope
        polyscope::init();
        meshV = simulation.currentV;
        meshF = simulation.faces;
        psMesh = polyscope::registerSurfaceMesh("sail", meshV, meshF);
        psMesh->setSmoothShade(true);
        psMesh->setSurfaceColor(glm::vec3(252. / 255., 247. / 255., 216. / 255.));
        psMesh->setEdgeWidth(0.0);

        // TODO: change place?
        Eigen::MatrixXd boatV;
        Eigen::MatrixXi boatF;
        igl::readOBJ("../../../Projects/PBD/data/sailboat/parts/boat.obj", boatV, boatF);
        polyscope::SurfaceMesh* boatMesh = polyscope::registerSurfaceMesh("boat", boatV, boatF);
        boatMesh->setSmoothShade(false);
        boatMesh->setSurfaceColor(glm::vec3(161. / 255., 77. / 255., 34. / 255.));
        boatMesh->setEdgeWidth(0.0);
        
        Eigen::MatrixXd waterV;
        Eigen::MatrixXi waterF;
        igl::readOBJ("../../../Projects/PBD/data/environment/water/water.obj", waterV, waterF);
        polyscope::SurfaceMesh* waterMesh = polyscope::registerSurfaceMesh("water", waterV, waterF);
        waterMesh->setSmoothShade(true);
        waterMesh->setTransparency(0.5);
        waterMesh->setMaterial("wax");
        waterMesh->setSurfaceColor(glm::vec3(15. / 255., 72. / 255., 110. / 255.));
        waterMesh->setEdgeWidth(0.0);

        Eigen::MatrixXd sandV;
        Eigen::MatrixXi sandF;
        igl::readOBJ("../../../Projects/PBD/data/environment/sand/sand.obj", sandV, sandF);
        polyscope::SurfaceMesh* sandMesh = polyscope::registerSurfaceMesh("sand", sandV, sandF);
        sandMesh->setSmoothShade(true);
        sandMesh->setMaterial("wax");
        sandMesh->setSurfaceColor(glm::vec3(253. / 255., 217. / 255., 141. / 255.));
        sandMesh->setEdgeWidth(0.0);

        Eigen::MatrixXd starfishV;
        Eigen::MatrixXi starfishF;
        igl::readOBJ("../../../Projects/PBD/data/environment/starfish/starfish.obj", starfishV, starfishF);
        polyscope::SurfaceMesh* starfishMesh = polyscope::registerSurfaceMesh("starfish", starfishV, starfishF);
        starfishMesh->setSmoothShade(true);
        starfishMesh->setMaterial("wax");
        starfishMesh->setSurfaceColor(glm::vec3(232. / 255., 124. / 255., 255. / 255.));
        starfishMesh->setEdgeWidth(0.0);

        Eigen::MatrixXd seaweedV;
        Eigen::MatrixXi seaweedF;
        igl::readOBJ("../../../Projects/PBD/data/environment/seaweed/seaweed_all.obj", seaweedV, seaweedF);
        polyscope::SurfaceMesh* seaweedMesh = polyscope::registerSurfaceMesh("seaweed", seaweedV, seaweedF);
        seaweedMesh->setSmoothShade(true);
        seaweedMesh->setMaterial("wax");
        seaweedMesh->setSurfaceColor(glm::vec3(46. / 255., 145. / 255., 31. / 255.));
        seaweedMesh->setEdgeWidth(0.0);

        Eigen::MatrixXd fishesV;
        Eigen::MatrixXi fishesF;
        igl::readOBJ("../../../Projects/PBD/data/environment/fishes/fishes_all.obj", fishesV, fishesF);
        polyscope::SurfaceMesh* fishesMesh = polyscope::registerSurfaceMesh("fishes", fishesV, fishesF);
        fishesMesh->setSmoothShade(true);
        fishesMesh->setMaterial("wax");
        fishesMesh->setSurfaceColor(glm::vec3(42. / 255., 135. / 255., 188. / 255.));
        fishesMesh->setEdgeWidth(0.0);

        Eigen::MatrixXd seashellV;
        Eigen::MatrixXi seashellF;
        igl::readOBJ("../../../Projects/PBD/data/environment/seashell/seashell_all.obj", seashellV, seashellF);
        polyscope::SurfaceMesh* seashellMesh = polyscope::registerSurfaceMesh("seashells", seashellV, seashellF);
        seashellMesh->setSmoothShade(true);
        seashellMesh->setMaterial("wax");
        seashellMesh->setSurfaceColor(glm::vec3(236. / 255., 226. / 255., 204. / 255.));
        seashellMesh->setEdgeWidth(0.0);

        glm::vec3 camPos = {-10.2436, 14.852, 5.19575};
        glm::vec3 camLookDir = {.481319, -0.8376, -0.258376};
        polyscope::view::lookAt(camPos, camLookDir);

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
        const double min_stiff = 0.0;
        const double max_stiff = 1.0;
        const double min_compliance = 0.0;
        const double max_compliance = 10.0;
        if (simulation.useXPBD) {
            ImGui::SliderScalar("Stretching compliance", ImGuiDataType_Double, &simulation.alpha_stretch, &min_compliance, &max_compliance, "%.3f");
            ImGui::SliderScalar("Bend compliance", ImGuiDataType_Double, &simulation.alpha_bend, &min_compliance, &max_compliance, "%.3f");
        } else {
            ImGui::SliderScalar("Stretching stiffness ", ImGuiDataType_Double, &simulation.k_stretch, &min_stiff, &max_stiff, "%.2f");
            ImGui::SliderScalar("Bend stiffness ", ImGuiDataType_Double, &simulation.k_bend, &min_stiff, &max_stiff, "%.2f");
        }

        ImGui::SliderScalar("Damping", ImGuiDataType_Double, &simulation.k_damping, &min_stiff, &max_stiff, "%.2f");

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
        ImGui::Checkbox("Fake Wind", &simulation.fakeWindActivated);
        ImGui::Checkbox("XPBD", &simulation.useXPBD);

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
