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
    double dt = 0.1;
    int solverPerRender = 10;
    MeshGenerator meshGen(1, Vec<double, 2>{0.0, DOMAIN_X_MAX}, numCells, MeshType::CARTESIAN);
    // MeshGenerator meshGen2(1, Vec<double, 2>{0.0, DOMAIN_X_MAX}, numCells, MeshType::CARTESIAN);
    Mesh<double, 1> mesh = meshGen.generateMesh();
    // Mesh<double, 1> mesh2 = meshGen2.generateMesh();
    Solver<double, 1> solver(mesh);
    // Solver<double, 1> solver2(mesh2);

    auto kernel = std::make_unique<WendlandFunction<double, 1>>(10 * DOMAIN_X_MAX/static_cast<double>(numCells));
    RBFs<double, 1> rbfs(
        numCells,
        numCells,
        mesh.cells.cellCentersPositions,
        std::move(kernel)
        );
    rbfs.computeLocalCellConnectivityMatrix(mesh);
    // rbfs.computeConnectivityMatrix(mesh);
    // std::cout << "Connectivity matrix computed." << std::endl;
    // float toleranceRBFExponent = -4.0;
    // uint maxIterationsRBF = 1000;
    // rbfs.computeCoefficientDerivatives(maxIterationsRBF*10, pow(10, toleranceRBFExponent));
    // std::cout << "Coefficient derivatives computed." << std::endl;
    
    auto doStep = [&]()
    {
        for (uint iter= 0; iter< solverPerRender; iter++)
{
            std::vector<double> plotY;
            // std::vector<double> plotY2;
            // std::vector<double> plotYRBFs;
            // plotX.reserve(mesh.cells.cellCentersPositions.size());
            plotY.reserve(mesh.cells.cellAverages.size());
            for (const auto& avg : mesh.cells.cellAverages)
            {
                plotY.push_back(avg[0]);
            }
            // rbfs.computeRBFCoefficients(plotY)
            // rbfs.computeRBFCoefficients(plotY, maxIterationsRBF, pow(10, toleranceRBFExponent));
            rbfs.computeRBFCoefficientsViaDerivatives(plotY);
            // std::cout << "RBF coefficients computed." << std::endl;
            solver.stepRBFBased(rbfs, dt);
            solver2.step(dt);
            // solver.step(dt);
            // std::cout << "Stepping"<<std::endl;
            // std::cout << "mesh average: "<< mesh.cells.cellAverages[10][0] <<std::endl;
        }
    };

    while(window.isOpen()){
        window.clear(sf::Color::Black);
         while (auto event = window.pollEvent()){
            ImGui::SFML::ProcessEvent(window, *event);
            if (event->is<sf::Event::Closed>())
            {
                window.close();
            }
             if (event->is<sf::Event::MouseButtonPressed>())
             {
                 doStep();
             }
        }
        ImGui::SFML::Update(window, clock.restart());
        ImGui::Begin("Linear advection");
        ImGui::SliderInt("Solver per render", &solverPerRender, 1, 1000);
        ImGui::SliderFloat("Tolerance exponent", &toleranceRBFExponent, -12.0, -1.0);
        ImGui::Text("Tolerance = %.5e", pow(10, toleranceRBFExponent));
        ImGui::SliderInt("Max iterations", reinterpret_cast<int*>(&maxIterationsRBF), 1, 10000);
        ImGui::Text("Iterations = %d", maxIterationsRBF);


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
            plotY2.reserve(mesh2.cells.cellAverages.size());
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
            for (const auto& avg : mesh2.cells.cellAverages)
            {
                plotY2.push_back(avg[0]);
            }
            ImPlot::SetupAxesLimits(0.0, DOMAIN_X_MAX, 0.0, 1.5);
            for (const auto& pos : mesh.cells.cellCentersPositions)
            {
                plotYRBFs.push_back(rbfs.value(pos));
            }
            // std::cout<<EnergyLog.data()<<std::endl;
            ImPlot::SetupAxes("x", "u");
            // std::cout<<"RBF plot [0] = "<<plotYRBFs[0]<<", [1] = "<<plotYRBFs[1]<<", [2] = "<<plotYRBFs[2]<<std::endl;
            // ImPlot::PlotLine("u", plotX.data(), plotY.data(), static_cast<int>(plotY.size()));
            ImPlot::PlotLine("u_RBFs", plotX.data(), plotYRBFs.data(), static_cast<int>(plotYRBFs.size()));
            ImPlot::PlotLine("u2", plotX.data(), plotY2.data(), static_cast<int>(plotY2.size()));
            ImPlot::EndPlot();
        }
        ImGui::End();
        ImGui::SFML::Render(window);
        window.display();
    }

    return 0;
}