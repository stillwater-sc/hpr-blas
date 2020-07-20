/*
// rot90.hpp : Rotate matrix 90 degrees.
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.
*/
//#include <boost/numeric/mtl/mtl.hpp>
//#include <hprblas>

#include <matpak/flipud.hpp>

//#define NOW
#ifdef NOW


namespace sw { namespace hprblas { namespace matpak {

template<typename Matrix>
Matrix rot90(Matrix&A) {
     return flipud(mtl::mat::trans(A));
}

}}} // namespace sw::hprblas::matpak
#endif
