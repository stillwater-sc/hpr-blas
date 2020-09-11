#pragma once
// rot90.hpp : Rotate matrix 90 degrees counter clockwise
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

#include <matpak/flipud.hpp>

namespace sw { namespace hprblas { namespace matpak {

template<typename Matrix>
Matrix rot90(const Matrix&A) {
     Matrix T(mtl::mat::trans(A));
     return flipud(T);
}
}}} 