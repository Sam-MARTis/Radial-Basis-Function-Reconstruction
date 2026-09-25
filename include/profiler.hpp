#pragma once

#pragma once
#include <unordered_map>
#include <vector>
#include  <string>
#include <chrono>
#include <iostream>
#include <algorithm>
#include <iomanip>

class Profiler
{

public:

    inline static std::unordered_map<std::string, double> logs;
    std::chrono::high_resolution_clock::time_point startTime;
    double logTimeInMS = 0;
    std::string logName;
    Profiler(const std::string& name);
    ~Profiler();
    static void printAndResetLogs();
};

inline Profiler::Profiler(const std::string& name)
{
    logName = name;
    startTime = std::chrono::high_resolution_clock::now();
}

inline Profiler::~Profiler()
{
    const std::chrono::high_resolution_clock::time_point destroyedTime = std::chrono::high_resolution_clock::now();
    logTimeInMS = std::chrono::duration<double, std::milli>(destroyedTime - startTime).count();
    logs[logName] += logTimeInMS;
}

inline void Profiler::printAndResetLogs()
{
    if (logs.empty())
        return;

    std::vector<std::pair<std::string, double>> entries(
        logs.begin(), logs.end());

    std::sort(entries.begin(), entries.end(),
        [](const auto& a, const auto& b)
        {
            return a.second > b.second;
        });

    double total = 0.0;
    for (const auto& [_, time] : entries)
        total += time;

    constexpr int nameWidth = 50;

    std::cout << "\n";
    std::cout << "==================== Profiling ====================\n";
    std::cout << std::left
              << std::setw(nameWidth) << "Section"
              << std::right
              << std::setw(12) << "Time (ms)"
              << std::setw(10) << "%\n";
    std::cout << "---------------------------------------------------\n";

    for (const auto& [name, time] : entries)
    {
        std::cout << std::left
                  << std::setw(nameWidth) << name
                  << std::right
                  << std::setw(12) << std::fixed << std::setprecision(3)
                  << time
                  << std::setw(9)
                  << std::setprecision(1)
                  << (100.0 * time / total)
                  << "%\n";
    }

    std::cout << "---------------------------------------------------\n";
    std::cout << std::left
              << std::setw(nameWidth) << "Total"
              << std::right
              << std::setw(12) << std::fixed << std::setprecision(3)
              << total << "\n\n";

    logs.clear();
}