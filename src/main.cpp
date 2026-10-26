#include <iostream>
#include <omp.h>
#include <SFML/Graphics.hpp>
#include <implot.h>

#include "imgui-SFML.h"

#include "mesh.hpp"
#include "fluid-solver.hpp"
#include "vec.hpp"
#include "RBF.hpp"
#include "validator.hpp"

// Conserved cell state -> primitive {rho, u, p} for plotting.
static Vec<double, 3> primitiveOf(const Vec<double, 3>& U)
{
    constexpr double gamma = 1.4;
    const double rho = std::max(U[0], EPSILON2);
    const double u = U[1] / rho;
    return Vec<double, 3>{rho, u, (gamma - 1.0) * (U[2] - 0.5 * rho * u * u)};
}


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
    double dt = 0.001;
    int solverPerRender = 20;
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

    // Conserved states {rho, rho*u, E}, consistent with EulerEquations / RoeFlux.
    Vec<double, 1+2> leftState{rhoL, rhoL*uL, EL};
    Vec<double, 1+2> rightState{rhoR, rhoR*uR, ER};

    Constraint<double, 1>::DirichletCondition(mesh, superDomain, rightState);
    Constraint<double, 1>::DirichletCondition(mesh, leftHalfDomain, leftState);


    // Mesh<double, 1> mesh2 = meshGen2.generateMesh();
    Solver<double, 1> solver(mesh);
    // Solver<double, 1> solver2(mesh2);

    solver.prescribeBoundaryConditions(0, BOUNDARY_CONSTRAINT::DIRICHLET, leftState);
    solver.prescribeBoundaryConditions(1, BOUNDARY_CONSTRAINT::DIRICHLET, rightState);
    std::cout<<"Boundary conditions prescribed."<<std::endl;

    RBF_KERNEL_KIND activeKernel = DEFAULT_RBF_KERNEL;
    float kernelSupportMultiplier = static_cast<float>(KERNEL_SUPPORT_MULTIPLIER);
    int limiterBracketDepth = static_cast<int>(LIMITER_NEIGHBOURHOOD_DEPTH);
    const auto kernelSupportRadius = [&]()
    {
        return static_cast<double>(kernelSupportMultiplier) * DEFAULT_NEIGHBOURHOOD_DEPTH *
               DOMAIN_X_MAX / static_cast<double>(numCells);
    };
    RBFs<double, 1> rbfs(
        numCells,
        numCells,
        mesh.cells.cellCentersPositions,
        makeRBFKernel<double, 1>(activeKernel, kernelSupportRadius())
        );
    rbfs.computeLocalCellConnectivityMatrix(mesh);
    std::cout<<"RBFs initialized with kernel "<<rbfKernelName(activeKernel)
             <<" (R = "<<rbfs.kernelSupportRadius()<<")."<<std::endl;

    const double interfacePosition = 0.5 * DOMAIN_X_MAX;
    riemann::ExactRiemannSolver<double> exactSolution(
        riemann::asPrimitive<double, 1>(leftState, riemann::VariableFormulation::CONSERVATIVE),
        riemann::asPrimitive<double, 1>(rightState, riemann::VariableFormulation::CONSERVATIVE),
        interfacePosition);
    std::cout << "Exact Riemann solution: p* = " << exactSolution.starPressure()
              << ", u* = " << exactSolution.starVelocity()
              << ", shock speed = " << exactSolution.tailSpeed();
    if (!exactSolution.converged()) std::cerr << " (NOT CONVERGED)";
    std::cout << std::endl;
    const riemann::SolutionValidator<double> validator(exactSolution, riemann::VariableFormulation::CONSERVATIVE);

    std::string name1 = "Euler Equations 1D ";
    std::string name2 = "RBF - Second Order";
    std::string name3 = "RBF - Direct";
    std::string name4 = "RBF - Second Order and Direct Blend";
    float simTime = 0;

    bool showExactSolution = true;
    bool showStatistics = true;

    // Drafted while only rho is plotted; other fields can be added here later.
    const riemann::FieldSpec densityField{riemann::FieldVariable::DENSITY, "rho", "kg/m^3", true};
    const std::vector<riemann::FieldSpec> validationFields{densityField};

    riemann::VolumetricData<double> validation;
    double lastValidatedTime = -1.0;
    bool haveStatistics = false;

    // Kernel changes reset the run.
    const auto resetSolution = [&]()
    {
        Constraint<double, 1>::DirichletCondition(mesh, superDomain, rightState);
        Constraint<double, 1>::DirichletCondition(mesh, leftHalfDomain, leftState);
        simTime = 0;
        lastValidatedTime = -1.0;
        haveStatistics = false;
    };
    const auto rebuildKernel = [&]()
    {
        rbfs.setKernel(makeRBFKernel<double, 1>(activeKernel, kernelSupportRadius()));
        rbfs.computeLocalCellConnectivityMatrix(mesh);
        resetSolution();
        std::cout << "Kernel switched to " << rbfKernelName(activeKernel)
                  << " (R = " << rbfs.kernelSupportRadius() << ")." << std::endl;
    };

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
        ImGui::Begin((name1 + name2).c_str());
        ImGui::SliderInt("Solver per render", &solverPerRender, 1, 1000);
        // ImGui::SliderFloat("Tolerance exponent", &toleranceRBFExponent, -12.0, -1.0);
        // ImGui::Text("Tolerance = %.5e", pow(10, toleranceRBFExponent));
        // ImGui::SliderInt("Max iterations", reinterpret_cast<int*>(&maxIterationsRBF), 1, 10000);
        // ImGui::Text("Iterations = %d", maxIterationsRBF);
        ImGui::Text("Simulation Time = %.2f", simTime);

        const char* kernelItems[] = {
            rbfKernelName(RBF_KERNEL_KIND::WENDLAND_C0),
            rbfKernelName(RBF_KERNEL_KIND::WENDLAND_C2),
            rbfKernelName(RBF_KERNEL_KIND::WENDLAND_C4),
            rbfKernelName(RBF_KERNEL_KIND::WENDLAND_C6),
            rbfKernelName(RBF_KERNEL_KIND::WU_C2),
            rbfKernelName(RBF_KERNEL_KIND::BUMP_CINF),
            rbfKernelName(RBF_KERNEL_KIND::WENDLAND_D4_C4),
            rbfKernelName(RBF_KERNEL_KIND::WU_C5),
            rbfKernelName(RBF_KERNEL_KIND::COMPACT_EXP),
            rbfKernelName(RBF_KERNEL_KIND::COMPACT_RATIONAL),
        };
        int kernelIndex = static_cast<int>(activeKernel);
        bool kernelJustChanged = false;
        if (ImGui::Combo("RBF kernel", &kernelIndex, kernelItems, IM_ARRAYSIZE(kernelItems)))
        {
            activeKernel = static_cast<RBF_KERNEL_KIND>(kernelIndex);
            rebuildKernel();
            kernelJustChanged = true;
        }
        if (ImGui::SliderFloat("Kernel support / (depth*h)", &kernelSupportMultiplier, 0.1f, 3.0f, "%.2f"))
        {
            rebuildKernel();
            kernelJustChanged = true;
        }
        const char* limiterItems[] = {
            limiterName(LIMITER_KIND::NO_LIMITER),
            limiterName(LIMITER_KIND::VENKATAKRISHNAN),
            limiterName(LIMITER_KIND::VAN_LEER),
            limiterName(LIMITER_KIND::MINMOD),
            limiterName(LIMITER_KIND::SUPERBEE),
            limiterName(LIMITER_KIND::MC),
            limiterName(LIMITER_KIND::VAN_ALBADA),
            limiterName(LIMITER_KIND::BARTH_JESPERSEN),
            limiterName(LIMITER_KIND::KOREN),
        };
        int limiterIndex = static_cast<int>(solver.limiter());
        if (ImGui::Combo("Limiter", &limiterIndex, limiterItems, IM_ARRAYSIZE(limiterItems)))
        {
            solver.setLimiter(static_cast<LIMITER_KIND>(limiterIndex));
            resetSolution();
            kernelJustChanged = true;
        }
        if (ImGui::SliderInt("Limiter bracket depth", &limiterBracketDepth, 1,
                             static_cast<int>(DEFAULT_NEIGHBOURHOOD_DEPTH)))
        {
            Mesh<double, 1>::reinitializeNeighbourhood(mesh, DEFAULT_NEIGHBOURHOOD_DEPTH,
                                                       static_cast<uint>(limiterBracketDepth));
            rbfs.computeLocalCellConnectivityMatrix(mesh);
            resetSolution();
            kernelJustChanged = true;
        }
        ImGui::Text("kernel: %s   R = %.4f (%.2f h)",
                    rbfKernelName(activeKernel),
                    rbfs.kernelSupportRadius(),
                    rbfs.kernelSupportRadius() / (DOMAIN_X_MAX / static_cast<double>(numCells)));
        ImGui::Text("limiter: %s   bracket: %d hop(s)   RBF stencil: %u hops",
                    limiterName(solver.limiter()), limiterBracketDepth,
                    DEFAULT_NEIGHBOURHOOD_DEPTH);
        if (ImGui::Button("Reset solution"))
        {
            resetSolution();
            kernelJustChanged = true;
        }
        ImGui::SameLine();
        if (ImGui::Button("Low-oscillation preset"))
        {
            activeKernel = RBF_KERNEL_KIND::WENDLAND_C4;
            kernelSupportMultiplier = 1.9f;
            limiterBracketDepth = 1;
            solver.setLimiter(LIMITER_KIND::MINMOD);
            Mesh<double, 1>::reinitializeNeighbourhood(mesh, DEFAULT_NEIGHBOURHOOD_DEPTH, 1);
            rebuildKernel();
            kernelJustChanged = true;
        }

        if (!kernelJustChanged && simTime < 3.0)
        {
            for (uint iter = 0; iter < solverPerRender; iter++)
            {
                solver.stepRBFBased(rbfs, dt);
                simTime += static_cast<float>(dt);
                // solver.stepPiecewise(dt);
            }
        }

        // ---- Exact-solution comparison controls ------------------------------
        ImGui::SeparatorText("Exact solution");
        ImGui::Checkbox("Show exact solution", &showExactSolution);
        ImGui::SameLine();
        ImGui::Checkbox("Show statistics", &showStatistics);

        // The validation walks every cell, so refresh it only when the solution changed.
        if (showExactSolution && (!haveStatistics ||
                                  std::abs(static_cast<double>(simTime) - lastValidatedTime) > 1e-9))
        {
            validation = validator.validate(mesh, static_cast<double>(simTime), validationFields);
            lastValidatedTime = static_cast<double>(simTime);
            haveStatistics = true;
        }

        if (showStatistics && haveStatistics)
        {
            const riemann::ConservationStats<double>& cons = validation.conservation;
            ImGui::SeparatorText("Error statistics (vs exact solution)");
            ImGui::Text("exact: p* = %.6g   u* = %.6g   shock speed = %.6g",
                        validation.exactInfo.pStar, validation.exactInfo.uStar, validation.exactInfo.shockSpeed);
            ImGui::Text("conservation drift:  mass %.3e (rel %.2e)   momentum %.3e (rel %.2e)   energy %.3e (rel %.2e)",
                        cons.massError, cons.massErrorRel,
                        cons.momentumError, cons.momentumErrorRel,
                        cons.energyError, cons.energyErrorRel);

            const ImGuiTableFlags tableFlags = ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg |
                                               ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_ScrollY;
            if (ImGui::BeginTable("validation_stats", 6, tableFlags, ImVec2(0.0f, 175.0f)))
            {
                ImGui::TableSetupScrollFreeze(0, 1);
                ImGui::TableSetupColumn("field");
                ImGui::TableSetupColumn("L1");
                ImGui::TableSetupColumn("L2");
                ImGui::TableSetupColumn("Linf");
                ImGui::TableSetupColumn("rel L1");
                ImGui::TableSetupColumn("integral error");
                ImGui::TableHeadersRow();
                for (const riemann::FieldStats<double>& field : validation.fields)
                {
                    ImGui::TableNextRow();
                    ImGui::TableNextColumn();
                    ImGui::TextUnformatted(field.name.c_str());
                    ImGui::TableNextColumn();
                    ImGui::Text("%.4e", field.l1);
                    ImGui::TableNextColumn();
                    ImGui::Text("%.4e", field.l2);
                    ImGui::TableNextColumn();
                    ImGui::Text("%.4e", field.lInf);
                    ImGui::TableNextColumn();
                    ImGui::Text("%.4e", field.relL1);
                    ImGui::TableNextColumn();
                    ImGui::Text("%.4e", field.errorIntegral);
                }
                ImGui::EndTable();
            }
        }

        // std::cout<<"Solver steps completed."<<std::endl;
        // ImPlot::SetNextAxesToFit();

        // rbfs.computeRBFCoefficients(mesh.cells.cellAverages);



        // doStep();

        if (ImPlot::BeginPlot("Euler Equation", ImVec2(-1, 500)))
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

            ImPlotSpec simulatedSpec;
            simulatedSpec.Marker = ImPlotMarker_Circle;
            simulatedSpec.MarkerSize = 2.5f;
            ImPlotSpec exactSpec;
            exactSpec.LineWeight = 2.0f;

            std::vector<double> exactY;
            if (showExactSolution) exactY.resize(plotX.size(), 0.0);

            ImPlot::SetupAxes("x", densityField.name);

            plotY.clear();
            for (const Vec<double, 3>& U : mesh.cells.cellAverages)
            {
                plotY.push_back(riemann::extractField<double, 1>(
                    densityField.variable, primitiveOf(U), U));
            }
            if (showExactSolution)
            {
                for (size_t cell = 0; cell < plotX.size(); ++cell)
                {
                    exactY[cell] = riemann::extractField<double, 1>(
                        densityField.variable,
                        exactSolution.samplePrimitive(plotX[cell], static_cast<double>(simTime)),
                        exactSolution.sampleConserved(plotX[cell], static_cast<double>(simTime)));
                }
            }

            ImPlot::PlotLine("rho simulated", plotX.data(), plotY.data(), static_cast<int>(plotY.size()), simulatedSpec);
            if (showExactSolution)
            {
                ImPlot::PlotLine("rho exact", plotX.data(), exactY.data(), static_cast<int>(exactY.size()), exactSpec);
            }
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