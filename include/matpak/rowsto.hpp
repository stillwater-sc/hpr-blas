#pragma once
// rowsto.hpp : Generates a random m x n row stochastic matrix
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.
#include <cstdint>
#include <random>

namespace sw { namespace hprblas { namespace matpak {

//
// size_t 64-bit unsigned
//
template<typename Matrix>
Matrix rowsto(size_t m, size_t n) {
    typedef typename Matrix::value_type value_type;

    // Use random_device to generate a seed for Mersenne twister engine.
    std::random_device rd{};
    // Use Mersenne twister engine to generate pseudo-random numbers.
    std::mt19937 engine{ rd() };
    // "Filter" MT engine's output to generate pseudo-random double values,
    // **uniformly distributed** on the closed interval [lowerbound, upperbound].
    // (Note that the range is [inclusive, inclusive].)
    constexpr double lowerbound = 0;
    constexpr double upperbound =  10;

    std::uniform_real_distribution<double> dist{ lowerbound, upperbound };
    // Pattern to generate pseudo-random number.
    // double rnd_value = dist(engine);

    Matrix A(m,n);
    int i ,j;

    for(i = 0; i < m; ++i) {
        for( j = 0;  j < n; ++j) {
            A[i][j] = value_type(dist(engine));
        }
        value_type x = mtl::sum(A[i][mtl::iall]);
        for( j = 0;  j < n; ++j) {
             A[i][j] = A[i][j]/x;
        }
    }
    return A;
}
}}} // namespace sw::hprblas::matpak
