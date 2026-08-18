#include <iostream>
#include <SFML/Graphics.hpp>
#include <implot.h>

#include "imgui-SFML.h"

#include "mesh.hpp"
#include "fluid-solver.hpp"
#include "vec.hpp"
#include "RBF.hpp"


int main()
{
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


    uint numCells = 200;
    double dt = 0.01;
    MeshGenerator meshGen(1, Vec<double, 2>{0.0, DOMAIN_X_MAX}, numCells, MeshType::CARTESIAN);
    Mesh<double, 1> mesh = meshGen.generateMesh();
    Solver<double, 1> solver(mesh.edges, mesh.cells, 0.1);

    auto kernel = std::make_unique<WendlandFunction<double, 1>>(3 * DOMAIN_X_MAX/static_cast<double>(numCells));
    RBFs<double, 1> rbfs(
        numCells,
        numCells,
        mesh.cells.cellCentersPositions,
        std::move(kernel)
        );

    rbfs.computeConnectivityMatrix(mesh);
    int solverPerRender = 30;


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
        ImGui::Begin("Linear advection");
        ImGui::SliderInt("Solver per render", &solverPerRender, 1, 1000);
        // ImPlot::SetNextAxesToFit();

        // rbfs.computeRBFCoefficients(mesh.cells.cellAverages);
        if (ImPlot::BeginPlot("Linear advection"))
        {
            std::vector<double> plotX;
            std::vector<double> plotY;
            std::vector<double> plotYRBFs;
            plotX.reserve(mesh.cells.cellCentersPositions.size());
            plotY.reserve(mesh.cells.cellAverages.size());
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
            rbfs.computeRBFCoefficients(plotY);
            for (const auto& pos : mesh.cells.cellCentersPositions)
            {
                plotYRBFs.push_back(rbfs.value(pos));
            }
            ImPlot::SetupAxesLimits(0.0, DOMAIN_X_MAX, 0.0, 1.5);

            // std::cout<<EnergyLog.data()<<std::endl;
            ImPlot::SetupAxes("x", "u");
            std::cout<<"RBF plot [0] = "<<plotYRBFs[0]<<", [1] = "<<plotYRBFs[1]<<", [2] = "<<plotYRBFs[2]<<std::endl;
            ImPlot::PlotLine("u", plotX.data(), plotY.data(), static_cast<int>(plotY.size()));
            ImPlot::PlotLine("u_RBFs", plotX.data(), plotYRBFs.data(), static_cast<int>(plotYRBFs.size()));
            ImPlot::EndPlot();
        }
        for (uint iter= 0; iter< solverPerRender; iter++)
        {
            solver.step(dt);
        }
        ImGui::End();
        ImGui::SFML::Render(window);
        window.display();
    }

    return 0;
}