#pragma once

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "constants.hpp"
#include "vec.hpp"

namespace sod
{

enum class Field
{
    DENSITY,
    VELOCITY,
    PRESSURE,
    MACH,
};

struct Conditions
{
    double gamma = 1.4;
    double time = 0.2;
    double domainLength = 1.0;
    double interfaceX = 0.5;
    double leftRho = 1.0;
    double leftU = 0.0;
    double leftP = 1.0;
    double rightRho = 0.125;
    double rightU = 0.0;
    double rightP = 0.1;
    double shockX = 0.850431;
    double contactX = 0.685490;
    std::string source;
    std::string sourceUrl;
    std::string reference;
};

inline std::string trim(const std::string& value)
{
    const auto begin = value.find_first_not_of(" \t\r\n");
    if (begin == std::string::npos) return {};
    const auto end = value.find_last_not_of(" \t\r\n");
    return value.substr(begin, end - begin + 1);
}

inline std::vector<std::string> split(const std::string& value, const char delimiter)
{
    std::vector<std::string> parts;
    std::stringstream stream(value);
    std::string part;
    while (std::getline(stream, part, delimiter)) parts.push_back(trim(part));
    return parts;
}

inline double parseNumber(const std::string& key, const std::string& value)
{
    try
    {
        std::size_t consumed = 0;
        const double result = std::stod(value, &consumed);
        if (consumed != value.size()) throw std::invalid_argument("trailing characters");
        return result;
    }
    catch (const std::exception&)
    {
        throw std::runtime_error("invalid numeric value for '" + key + "': " + value);
    }
}

class Data
{
public:
    static Data load(const std::filesystem::path& path)
    {
        std::ifstream input(path);
        if (!input) throw std::runtime_error("cannot open Sod data: " + path.string());

        Data data;
        std::unordered_map<std::string, std::string> metadata;
        std::vector<std::string> header;
        std::vector<std::vector<double>> columns;
        std::string line;
        bool headerRead = false;

        while (std::getline(input, line))
        {
            line = trim(line);
            if (line.empty()) continue;
            if (line.front() == '#')
            {
                const std::string entry = trim(line.substr(1));
                const auto equals = entry.find('=');
                if (equals != std::string::npos)
                {
                    metadata[trim(entry.substr(0, equals))] = trim(entry.substr(equals + 1));
                }
                continue;
            }
            if (!headerRead)
            {
                header = split(line, ',');
                columns.resize(header.size());
                headerRead = true;
                continue;
            }
            const std::vector<std::string> fields = split(line, ',');
            if (fields.size() != header.size())
            {
                throw std::runtime_error("column count mismatch in " + path.string());
            }
            for (std::size_t column = 0; column < fields.size(); ++column)
            {
                columns[column].push_back(parseNumber(header[column], fields[column]));
            }
        }

        if (!headerRead) throw std::runtime_error("missing CSV header in " + path.string());

        const auto number = [&metadata](const std::string& key)
        {
            const auto found = metadata.find(key);
            if (found == metadata.end()) throw std::runtime_error("missing metadata key: " + key);
            return parseNumber(key, found->second);
        };
        const auto text = [&metadata](const std::string& key)
        {
            const auto found = metadata.find(key);
            if (found == metadata.end()) throw std::runtime_error("missing metadata key: " + key);
            return found->second;
        };

        Conditions& conditions = data.conditions_;
        conditions.gamma = number("gamma");
        conditions.time = number("time_s");
        conditions.domainLength = number("domain_length_m");
        conditions.interfaceX = number("interface_x_m");
        conditions.leftRho = number("left_rho");
        conditions.leftU = number("left_u");
        conditions.leftP = number("left_p");
        conditions.rightRho = number("right_rho");
        conditions.rightU = number("right_u");
        conditions.rightP = number("right_p");
        conditions.shockX = number("shock_x_m");
        conditions.contactX = number("contact_x_m");
        conditions.source = text("source");
        conditions.sourceUrl = text("source_url");
        conditions.reference = text("reference");

        const auto column = [&header](const std::string& name)
        {
            const auto found = std::find(header.begin(), header.end(), name);
            if (found == header.end()) throw std::runtime_error("missing column: " + name);
            return static_cast<std::size_t>(found - header.begin());
        };

        const std::vector<double>& x = columns[column("x")];
        data.x_ = x;
        data.rho_ = columns[column("rho")];
        data.u_ = columns[column("u")];
        data.p_ = columns[column("p")];
        data.mach_ = columns[column("mach")];

        for (std::size_t i = 1; i < data.x_.size(); ++i)
        {
            if (!(data.x_[i] > data.x_[i - 1])) throw std::runtime_error("x column is not increasing");
        }
        if (data.x_.empty()) throw std::runtime_error("empty Sod data");

        return data;
    }

    const Conditions& conditions() const { return conditions_; }
    const std::vector<double>& x() const { return x_; }
    const std::vector<double>& rho() const { return rho_; }
    const std::vector<double>& u() const { return u_; }
    const std::vector<double>& p() const { return p_; }
    const std::vector<double>& mach() const { return mach_; }
    std::size_t size() const { return x_.size(); }

    double value(const Field field, const std::size_t index) const
    {
        switch (field)
        {
            case Field::DENSITY: return rho_.at(index);
            case Field::VELOCITY: return u_.at(index);
            case Field::PRESSURE: return p_.at(index);
            case Field::MACH: return mach_.at(index);
        }
        return 0.0;
    }

    double soundSpeedLeft() const
    {
        return std::sqrt(conditions_.gamma * conditions_.leftP / conditions_.leftRho);
    }

    double soundSpeedRight() const
    {
        return std::sqrt(conditions_.gamma * conditions_.rightP / conditions_.rightRho);
    }

    double shockSpeed() const
    {
        return (conditions_.shockX - conditions_.interfaceX) / conditions_.time;
    }

    double maxWaveSpeed() const
    {
        return std::max({soundSpeedLeft(), soundSpeedRight(), std::abs(shockSpeed())});
    }

    Vec<double, 3> leftConserved() const
    {
        return {conditions_.leftRho,
                conditions_.leftRho * conditions_.leftU,
                conditions_.leftP / (conditions_.gamma - 1.0)};
    }

    Vec<double, 3> rightConserved() const
    {
        return {conditions_.rightRho,
                conditions_.rightRho * conditions_.rightU,
                conditions_.rightP / (conditions_.gamma - 1.0)};
    }

private:
    Conditions conditions_;
    std::vector<double> x_;
    std::vector<double> rho_;
    std::vector<double> u_;
    std::vector<double> p_;
    std::vector<double> mach_;
};

inline std::filesystem::path defaultDataPath()
{
#ifdef SOD_DATA_PATH
    const std::filesystem::path compiled = SOD_DATA_PATH;
    if (std::filesystem::exists(compiled)) return compiled;
#endif
    for (const char* candidate : {"data/sod_shock_tube.csv",
                                  "../data/sod_shock_tube.csv",
                                  "../../data/sod_shock_tube.csv"})
    {
        if (std::filesystem::exists(candidate)) return candidate;
    }
    throw std::runtime_error("Sod CSV not found; pass its path explicitly");
}

}
