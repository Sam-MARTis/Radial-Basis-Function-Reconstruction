#pragma once
#include "SFML/Graphics.hpp"
#include <random>


template <typename T>
constexpr T absVal(T value) {
    return (value < 0) ? -value : value;
}

inline float randf(const float min, const float max)
{
    return min + (static_cast<float>(rand())/(static_cast<float>(RAND_MAX))) * (max-min);
}

struct HSL
{
    float h; // [0, 360)
    float s; // [0, 1]
    float l; // [0, 1]
};

inline sf::Color HSLToRGB(const HSL& hsl)
{
    const float h = hsl.h / 360.0f;
    const float s = hsl.s;
    const float l = hsl.l;

    if (s == 0.0f)
    {
        const auto v = static_cast<std::uint8_t>(l * 255.0f);
        return sf::Color(v, v, v, 255);
    }

    const float q = (l < 0.5f)
        ? l * (1.0f + s)
        : l + s - l * s;

    const float p = 2.0f * l - q;

    auto hueToRGB = [](float p, float q, float t)
    {
        if (t < 0.0f) t += 1.0f;
        if (t > 1.0f) t -= 1.0f;

        if (t < 1.0f / 6.0f)
            return p + (q - p) * 6.0f * t;

        if (t < 1.0f / 2.0f)
            return q;

        if (t < 2.0f / 3.0f)
            return p + (q - p) * (2.0f / 3.0f - t) * 6.0f;

        return p;
    };

    const float r = hueToRGB(p, q, h + 1.0f / 3.0f);
    const float g = hueToRGB(p, q, h);
    const float b = hueToRGB(p, q, h - 1.0f / 3.0f);

    return sf::Color(
        static_cast<std::uint8_t>(r * 255.0f),
        static_cast<std::uint8_t>(g * 255.0f),
        static_cast<std::uint8_t>(b * 255.0f),
        255
    );
}