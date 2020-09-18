#pragma once
// gt.hpp : A > B (determines if A > B element-wise)
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.
#include<matpak/relop.hpp>

namespace sw { namespace hprblas { namespace matpak {

template<typename Matrix>
Matrix operator>(const Matrix&A , const Matrix&B){
	return relop(A,RELOP_GT,B);
}

template<typename Matrix>
Matrix gt(const Matrix&A, const Matrix&B) {
	return relop(A,RELOP_GT,B);
}
}}}