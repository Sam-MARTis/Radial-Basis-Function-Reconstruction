#pragma once
#include "function.hpp"






template<typename T, uint N>
class RBFs{
    uint numRBFs = 0;
    Function<T, N> kernelFunction;
    std::vector<Vec<T, N>> kernelCenters;

};


// class MultiQuadratic