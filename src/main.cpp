#include <iostream>
#include "mesh.hpp"
#include "fluid-solver.hpp"
#include "vec.hpp"
#include <SFML/Graphics.hpp>
#include "imgui-SFML.h"
#include <implot.h>

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


    uint numCells = 500;
    double dt = 0.01;
    MeshGenerator meshGen(1, Vec<double, 2>{0.0, DOMAIN_X_MAX}, numCells, MeshType::CARTESIAN);
    Mesh<double, 1> mesh = meshGen.generateMesh();
    Solver<double, 1> solver(mesh.edges, mesh.cells, 0.1);



    int solverPerRender = 10;


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
        if (ImPlot::BeginPlot("Linear advection"))
        {
            std::vector<double> plotX;
            std::vector<double> plotY;
            plotX.reserve(mesh.cells.cellCentersPositions.size());
            plotY.reserve(mesh.cells.cellAverages.size());
            assert(mesh.cells.cellCentersPositions.size() == mesh.cells.cellAverages.size());

            for (const auto& pos : mesh.cells.cellCentersPositions)
            {
                plotX.push_back(pos[0]);
            }
            for (const auto& avg : mesh.cells.cellAverages)
            {
                plotY.push_back(avg[0]);
            }
            ImPlot::SetupAxesLimits(0.0, DOMAIN_X_MAX, 0.0, 1.5);

            // std::cout<<EnergyLog.data()<<std::endl;
            ImPlot::SetupAxes("x", "u");

            ImPlot::PlotLine("u", plotX.data(), plotY.data(), static_cast<int>(plotY.size()));
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