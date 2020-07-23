#prama once

// hankel.hpp : Generate Hankel matrix
// CALL hankel(c,r) where c = column & r = row
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

#include <matpak/toeplitz.hpp>
#include <matpak/rot90.hpp>

namespace sw { namespace hprblas { namespace matpak {

template<typename Vector>
mtl::dense2D<typename Vector::value_type> hankel(const Vector&c, const Vector&r) {
    typedef typename Vector::value_type value_type;

    using Matrix = mtl::dense2D<value_type>;

    auto A = toeplitz(c,r);

    // std::cout << A << std::endl;
    // auto B = rot90(A);


    return A; // return Hankel matrix
}
}}} // namespace sw::hprblas::matpak
