/*
// toeplitz.hpp : Generate Toeplitz matrix
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.
*/
#include <vector>

namespace sw { namespace hprblas { namespace matpak {

template<typename Vector>
mtl::dense2D<typename Vector::value_type> toeplitz(const Vector&c, const Vector&r) {
    typedef typename Vector::value_type value_type;

    using Matrix = mtl::dense2D<value_type>;

    Matrix A(4,4);
    return A;
}
}}} // namespace sw::hprblas::matpak
