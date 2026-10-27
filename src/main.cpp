#include <filesystem>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <omp.h>
#include <SFML/Graphics.hpp>
#include <implot.h>

#include "imgui-SFML.h"

#include "sod_shock_tube.hpp"
#include "sod_simulation.hpp"
#include "sod_validation.hpp"

int main(int argc, char** argv)
{
    omp_set_num_threads(std::max(1, omp_get_max_threads() - 1));

    std::filesystem::path dataPath;
    sod::Data data;
    try
    {
        dataPath = argc > 1 ? std::filesystem::path(argv[1]) : sod::defaultDataPath();
        data = sod::Data::load(dataPath);
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    sod::Options options;
    options.cells = 200;

    sf::ContextSettings settings;
    settings.antiAliasingLevel = 16;
    sf::RenderWindow window(sf::VideoMode({static_cast<unsigned int>(DEFAULT_SCREEN_WIDTH),
                                           static_cast<unsigned int>(DEFAULT_SCREEN_HEIGHT)}),
                            "Sod shock tube validation");
    if (!ImGui::SFML::Init(window)) return 1;
    ImPlot::CreateContext();
    window.setFramerateLimit(60);

    std::unique_ptr<sod::Simulation> simulation;
    sod::Validation validation;
    std::string status = "ready";
    std::vector<double> cellX;
    std::vector<double> density;
    std::vector<double> velocity;
    std::vector<double> pressure;

    const auto rebuild = [&]()
    {
        try
        {
            auto candidate = std::make_unique<sod::Simulation>(data, options);
            candidate->runToEnd();
            auto result = sod::Validator::validate(data, *candidate);

            cellX.resize(candidate->mesh().numCells);
            density.resize(candidate->mesh().numCells);
            velocity.resize(candidate->mesh().numCells);
            pressure.resize(candidate->mesh().numCells);
            for (uint cell = 0; cell < candidate->mesh().numCells; ++cell)
            {
                cellX[cell] = candidate->mesh().cells.cellCentersPositions[cell][0];
                density[cell] = candidate->value(sod::Field::DENSITY, cell);
                velocity[cell] = candidate->value(sod::Field::VELOCITY, cell);
                pressure[cell] = candidate->value(sod::Field::PRESSURE, cell);
            }

            simulation = std::move(candidate);
            validation = std::move(result);
            status = "ok";
        }
        catch (const std::exception& error)
        {
            simulation.reset();
            validation = {};
            status = error.what();
        }
    };
    rebuild();

    bool showReference = true;
    sf::Clock clock;

    while (window.isOpen())
    {
        while (auto event = window.pollEvent())
        {
            ImGui::SFML::ProcessEvent(window, *event);
            if (event->is<sf::Event::Closed>()) window.close();
        }
        ImGui::SFML::Update(window, clock.restart());
        ImGui::SetNextWindowPos(ImVec2(10.0f, 10.0f), ImGuiCond_Always);
        ImGui::SetNextWindowSize(ImVec2(static_cast<float>(DEFAULT_SCREEN_WIDTH) - 20.0f,
                                        static_cast<float>(DEFAULT_SCREEN_HEIGHT) - 30.0f),
                                 ImGuiCond_Always);
        ImGui::Begin("Sod shock tube validation");

        ImGui::Text("data: %s", dataPath.string().c_str());
        ImGui::Text("source: %s", data.conditions().source.c_str());
        ImGui::TextWrapped("%s", data.conditions().sourceUrl.c_str());
        ImGui::Text("reference: %s", data.conditions().reference.c_str());
        ImGui::Text("status: %s", status.c_str());
        ImGui::Separator();

        ImGui::Text("gamma = %.3f   t = %.4f   x0 = %.3f   cells = %u   dt = %.3e s",
                    data.conditions().gamma, data.conditions().time, data.conditions().interfaceX,
                    options.cells, simulation ? simulation->dt() : 0.0);
        ImGui::Text("left  rho/u/p = %.4g / %.4g / %.4g     right rho/u/p = %.4g / %.4g / %.4g",
                    data.conditions().leftRho, data.conditions().leftU, data.conditions().leftP,
                    data.conditions().rightRho, data.conditions().rightU, data.conditions().rightP);
        if (validation.valid)
        {
            ImGui::Text("min rho = %.4f   min p = %.6f   max Mach = %.4f   mass drift = %.3e",
                        validation.minimumDensity, validation.minimumPressure, validation.maximumMach,
                        validation.massErrorRel);
        }

        int cells = static_cast<int>(options.cells);
        if (ImGui::SliderInt("Cells", &cells, 50, 800))
        {
            if (ImGui::IsItemDeactivatedAfterEdit())
            {
                options.cells = static_cast<uint>(cells);
                rebuild();
            }
        }

        const char* kernels[] = {
            "Wendland C0", "Wendland C2", "Wendland C4", "Wendland C6", "Wu C2",
            "Bump C-inf", "Wendland D4C4", "Wu C5", "Compact exp", "Compact rational"};
        int kernel = static_cast<int>(options.kernel);
        if (ImGui::Combo("RBF kernel", &kernel, kernels, IM_ARRAYSIZE(kernels)))
        {
            options.kernel = static_cast<RBF_KERNEL_KIND>(kernel);
            rebuild();
        }

        float support = static_cast<float>(options.supportMultiplier);
        if (ImGui::SliderFloat("Support / (depth h)", &support, 0.1f, 3.0f, "%.2f"))
        {
            if (ImGui::IsItemDeactivatedAfterEdit())
            {
                options.supportMultiplier = static_cast<double>(support);
                rebuild();
            }
        }

        const char* limiters[] = {
            "none", "Venkatakrishnan", "van Leer", "minmod", "superbee",
            "MC", "van Albada", "Barth-Jespersen", "Koren"};
        int limiter = static_cast<int>(options.limiter);
        if (ImGui::Combo("Limiter", &limiter, limiters, IM_ARRAYSIZE(limiters)))
        {
            options.limiter = static_cast<LIMITER_KIND>(limiter);
            rebuild();
        }

        int bracket = static_cast<int>(options.limiterDepth);
        if (ImGui::SliderInt("Limiter bracket depth", &bracket, 1,
                             static_cast<int>(DEFAULT_NEIGHBOURHOOD_DEPTH)))
        {
            if (ImGui::IsItemDeactivatedAfterEdit())
            {
                options.limiterDepth = static_cast<uint>(bracket);
                rebuild();
            }
        }

        float cfl = static_cast<float>(options.cfl);
        if (ImGui::SliderFloat("CFL", &cfl, 0.05f, 0.5f, "%.2f"))
        {
            if (ImGui::IsItemDeactivatedAfterEdit())
            {
                options.cfl = static_cast<double>(cfl);
                rebuild();
            }
        }

        if (ImGui::Button("Recompute")) rebuild();
        ImGui::SameLine();
        ImGui::Checkbox("Reference data", &showReference);

        if (validation.valid)
        {
            ImGui::SeparatorText("Errors against the NASA reference profile");
            const ImGuiTableFlags flags = ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg |
                                          ImGuiTableFlags_SizingStretchProp;
            if (ImGui::BeginTable("errors", 6, flags))
            {
                ImGui::TableSetupColumn("field");
                ImGui::TableSetupColumn("L1");
                ImGui::TableSetupColumn("L2");
                ImGui::TableSetupColumn("Linf");
                ImGui::TableSetupColumn("rel L1");
                ImGui::TableSetupColumn("bias");
                ImGui::TableHeadersRow();
                for (const sod::FieldError& field : validation.fields)
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
                    ImGui::Text("%+.4e", field.bias);
                }
                ImGui::EndTable();
            }

            const auto field = [&](const char* name) -> const sod::FieldError&
            {
                return *validation.find(name);
            };

            if (ImGui::BeginTable("profiles", 3, ImGuiTableFlags_SizingStretchSame))
            {
                const sod::FieldError* fields[] = {&field("rho"), &field("u"), &field("p")};
                const std::vector<double>* series[] = {&density, &velocity, &pressure};
                const char* titles[] = {"Density", "Velocity", "Pressure"};
                const char* ylabels[] = {"rho", "u", "p"};

                for (int plot = 0; plot < 3; ++plot)
                {
                    ImGui::TableNextColumn();
                    if (ImPlot::BeginPlot(titles[plot], ImVec2(-1, 330)))
                    {
                        ImPlot::SetupAxes("x", ylabels[plot]);
                        ImPlotSpec line;
                        line.LineWeight = 1.6f;
                        ImPlot::PlotLine("simulation", cellX.data(), series[plot]->data(),
                                         static_cast<int>(cellX.size()), line);
                        if (showReference)
                        {
                            ImPlotSpec points;
                            points.Marker = ImPlotMarker_Circle;
                            points.MarkerSize = 4.0f;
                            ImPlot::PlotScatter("reference", fields[plot]->x.data(),
                                                fields[plot]->reference.data(),
                                                static_cast<int>(fields[plot]->x.size()), points);
                        }
                        ImPlot::EndPlot();
                    }
                }
                ImGui::EndTable();
            }

            ImGui::Text("conservation drift: mass %.3e (rel %.2e)   momentum %.3e (rel %.2e)   energy %.3e (rel %.2e)",
                        validation.massError, validation.massErrorRel,
                        validation.momentumError, validation.momentumErrorRel,
                        validation.energyError, validation.energyErrorRel);
        }
        else
        {
            ImGui::SeparatorText("Errors against the NASA reference profile");
            ImGui::TextUnformatted("No validation available.");
        }

        ImGui::End();
        ImGui::SFML::Render(window);
        window.display();
    }

    ImPlot::DestroyContext();
    ImGui::SFML::Shutdown();
    return 0;
}
