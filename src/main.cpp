#include <iostream>
#include <omp.h>
#include <SFML/Graphics.hpp>
#include <implot.h>

#include "imgui-SFML.h"

#include "mesh.hpp"
#include "fluid-solver.hpp"
#include "vec.hpp"
#include "RBF.hpp"


int main()
{

    std::cout<<"Max threads available: "<<omp_get_max_threads()<<std::endl;
    omp_set_num_threads(omp_get_max_threads()-1);


    sf::Clock clock;
    sf::ContextSettings contextSettings;
    contextSettings.antiAliasingLevel = 16;
    sf::RenderWindow window = sf::RenderWindow(sf::VideoMode({static_cast<unsigned int>(round(DEFAULT_SCREEN_WIDTH)), static_cast<unsigned int>(round(DEFAULT_SCREEN_HEIGHT))}),"FVM");
    bool imguiInitSuccess = ImGui::SFML::Init(window);
    if (!imguiInitSuccess)
        std::cerr << "SFML Init failed." << std::endl;
    ImPlot::CreateContext();
    window.setFramerateLimit(60);


    std::cout << "RenderWindow setup successful." << std::endl;


    uint numCells = 100;
    double dt = 0.01;
    int solverPerRender = 10;
    MeshGenerator meshGen(1, Vec<double, 2>{0.0, DOMAIN_X_MAX}, numCells, MeshType::CARTESIAN);
    // MeshGenerator meshGen2(1, Vec<double, 2>{0.0, DOMAIN_X_MAX}, numCells, MeshType::CARTESIAN);
    Mesh<double, 1> mesh = meshGen.generateMesh();
    mesh.fixEdgesConnectivityOrder();
    Vec<double, 1> center = (mesh.boundingBox.first + mesh.boundingBox.second)*0.5;
    Vec<double, 1> halfSpan = (mesh.boundingBox.second - mesh.boundingBox.first)*0.5;

    HyperCuboid<double, 1> leftBoundaryBoundingBox(center - halfSpan*1.1, center - halfSpan*0.9);
    HyperCuboid<double, 1> rightBoundaryBoundingBox(center + halfSpan*0.9, center + halfSpan*1.1);
    HyperCuboid<double, 1> superDomain(center - halfSpan*1.1, center + halfSpan*1.1);
    HyperCuboid<double, 1> leftHalfDomain(center - halfSpan*1.1, center);
    HyperCuboid<double, 1> rightHalfDomain(center, center + halfSpan*1.1);

    std::vector<std::pair<uint, std::unique_ptr<Domain<double, 1>>>> boundaryDomains;
    boundaryDomains.emplace_back(0, std::make_unique<HyperCuboid<double, 1>>(leftBoundaryBoundingBox));
    boundaryDomains.emplace_back(1, std::make_unique<HyperCuboid<double, 1>>(rightBoundaryBoundingBox));
    #define DO_STEP
    mesh.identifyBoundaryEdgesAndCells(boundaryDomains);
    const double rhoL = 1.0;
    const double rhoR = 0.125;
    const double pL = 1.0;
    const double pR = 0.1;
    const double uL = 0.05;
    const double uR = 0.0;
    const double EL = EulerEquations<double, 1>::getEFromPressureRhoVMagSq(pL, rhoL, uL*uL);
    const double ER = EulerEquations<double, 1>::getEFromPressureRhoVMagSq(pR, rhoR, uR*uR);

    Vec<double, 1+2> leftState{1.0, 0.0, 1.0};
    Vec<double, 1+2> rightState{0.125, 0.0, 0.1};

    Constraint<double, 1>::DirichletCondition(mesh, superDomain, rightState);
    Constraint<double, 1>::DirichletCondition(mesh, leftHalfDomain, leftState);


    // Mesh<double, 1> mesh2 = meshGen2.generateMesh();
    Solver<double, 1> solver(mesh);
    // Solver<double, 1> solver2(mesh2);

    solver.prescribeBoundaryConditions(0, BOUNDARY_CONSTRAINT::DIRICHLET, leftState);
    solver.prescribeBoundaryConditions(1, BOUNDARY_CONSTRAINT::DIRICHLET, rightState);
    std::cout<<"Boundary conditions prescribed."<<std::endl;

    auto kernel = std::make_unique<WendlandFunction<double, 1>>(10 * DOMAIN_X_MAX/static_cast<double>(numCells));
    RBFs<double, 1> rbfs(
        numCells,
        numCells,
        mesh.cells.cellCentersPositions,
        std::move(kernel)
        );
    rbfs.computeLocalCellConnectivityMatrix(mesh);
    std::cout<<"RBFs initialized."<<std::endl;

    while(window.isOpen()){
        window.clear(sf::Color::Black);
         while (auto event = window.pollEvent()){
            ImGui::SFML::ProcessEvent(window, *event);
            if (event->is<sf::Event::Closed>())
            {
                window.close();
            }
        }
        ImGui::SFML::Update(window, clock.restart());
        ImGui::Begin("Euler Equations");
        ImGui::SliderInt("Solver per render", &solverPerRender, 1, 1000);
        // ImGui::SliderFloat("Tolerance exponent", &toleranceRBFExponent, -12.0, -1.0);
        // ImGui::Text("Tolerance = %.5e", pow(10, toleranceRBFExponent));
        // ImGui::SliderInt("Max iterations", reinterpret_cast<int*>(&maxIterationsRBF), 1, 10000);
        // ImGui::Text("Iterations = %d", maxIterationsRBF);
        for (uint iter = 0; iter < solverPerRender; iter++)
        {
            // solver.stepRBFBased(rbfs, dt);
            solver.stepPiecewise(dt);
        }

        std::cout<<"Solver steps completed."<<std::endl;
        // ImPlot::SetNextAxesToFit();

        // rbfs.computeRBFCoefficients(mesh.cells.cellAverages);



        // doStep();

        if (ImPlot::BeginPlot("Linear advection"))
        {
            std::vector<double> plotX;
            std::vector<double> plotY;
            std::vector<double> plotY2;
            std::vector<double> plotYRBFs;
            plotX.reserve(mesh.cells.cellCentersPositions.size());
            plotY.reserve(mesh.cells.cellAverages.size());
            // plotY2.reserve(mesh2.cells.cellAverages.size());
            plotYRBFs.reserve(mesh.cells.cellAverages.size());
            assert(mesh.cells.cellCentersPositions.size() == mesh.cells.cellAverages.size());
            for (const auto& pos : mesh.cells.cellCentersPositions)
            {
                plotX.push_back(pos[0]);
            }
            for (const auto& avg : mesh.cells.cellAverages)
            {
                plotY.push_back(avg[0]);
            }
            // for (const auto& avg : mesh2.cells.cellAverages)
            // {
            //     plotY2.push_back(avg[0]);
            // }
            ImPlot::SetupAxesLimits(0.0, DOMAIN_X_MAX, 0.0, 1.5);
            // for (const auto& pos : mesh.cells.cellCentersPositions)
            // {
            //     plotYRBFs.push_back(rbfs.value(pos));
            // }
            // std::cout<<EnergyLog.data()<<std::endl;
            ImPlot::SetupAxes("x", "rho");
            // std::cout<<"RBF plot [0] = "<<plotYRBFs[0]<<", [1] = "<<plotYRBFs[1]<<", [2] = "<<plotYRBFs[2]<<std::endl;
            ImPlot::PlotLine("rho", plotX.data(), plotY.data(), static_cast<int>(plotY.size()));
            // ImPlot::PlotLine("u_RBFs", plotX.data(), plotYRBFs.data(), static_cast<int>(plotYRBFs.size()));
            // ImPlot::PlotLine("u2", plotX.data(), plotY2.data(), static_cast<int>(plotY2.size()));
            ImPlot::EndPlot();
        }
        ImGui::End();
        ImGui::SFML::Render(window);
        window.display();
    }

    return 0;
}

/* Reference:
 *    for (size_t i = 0; i < n; ++i)
    {

        const auto &pos = testMesh1D.cells.cellCentersPositions[i];
        const auto &U = testMesh1D.cells.cellAverages[i];
        plotX.push_back(pos[0]);
        const double r = U[0];
        const double u = (r != 0.0) ? (U[1] / r) : 0.0;
        const double E = U[2];
        const double p = (gamma - 1.0) * (E - 0.5 * r * u * u);
        rho.push_back(r);
        vel.push_back(u);
        pres.push_back(p);
    }

    if (ImGui::BeginTable("plots", 3, ImGuiTableFlags_SizingStretchProp))
    {
        ImGui::TableNextColumn();
        if (ImPlot::BeginPlot("Density", ImVec2(-1, 260)))
        {
            ImPlot::SetupAxes("x", "rho");
            ImPlot::SetupAxesLimits(0,1.0,0,1.5);
            ImPlot::PlotLine("rho", plotX.data(), rho.data(), static_cast<int>(n));
            ImPlot::EndPlot();
        }

        ImGui::TableNextColumn();
        if (ImPlot::BeginPlot("Velocity", ImVec2(-1, 260)))
        {
            ImPlot::SetupAxes("x", "u");
            ImPlot::SetupAxesLimits(0,1.0,0,1.25);
            ImPlot::PlotLine("u", plotX.data(), vel.data(), static_cast<int>(n));
            ImPlot::EndPlot();
        }

        ImGui::TableNextColumn();
        if (ImPlot::BeginPlot("Pressure", ImVec2(-1, 260)))
        {
            ImPlot::SetupAxes("x", "p");
            ImPlot::SetupAxesLimits(0,1.0,0,1.5);
            ImPlot::PlotLine("p", plotX.data(), pres.data(), static_cast<int>(n));
            ImPlot::EndPlot();
        }

        ImGui::EndTable();
    }

    for (uint iter= 0; iter< solverPerRender; iter++)
    {
        solver.step(dt);
    }
        // std::cout <<"Time step: " << dt << std::endl;
        if (iterationCount%100 == 0)
        {
            auto val = testMesh1D.cells.cellAverages[10][1];
            std::cout << "Iteration: " << iterationCount << ", cell 10 velocity: ";
            std::cout << val << std::endl;

        }
    ImGui::End();
    ImGui::SFML::Render(window);
        // testMeshRenderer.draw(window);
    window.display();
    }

    return 0;
}
*/