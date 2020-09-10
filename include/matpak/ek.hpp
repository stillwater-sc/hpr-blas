#pragma once
// ek.hpp : i-th standard basis vector of length n
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

namespace sw { namespace hprblas { namespace matpak {

template<typename Scalar>
mtl::vec::dense_vector<Scalar> ek(size_t i, size_t n ) {
    mtl::vec::dense_vector<Scalar> v(n, Scalar(0));
	v[i] = 1;
    return v;
}
}}}