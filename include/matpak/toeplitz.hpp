#pragma once
// toeplitz.hpp : Generate Toeplitz matrix
// CALL toeplitz(c,r) where c = column & r = row
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

namespace sw { namespace hprblas { namespace matpak {

template<typename Vector>
mtl::dense2D<typename Vector::value_type> toeplitz(const Vector&c, const Vector&r) {
    typedef typename Vector::value_type value_type;

    using Matrix = mtl::dense2D<value_type>;

    size_t m = num_rows(c);
	size_t n = num_rows(r);

	if (c[0]!=r[0]){
		std::cout <<
		"Warning: First element of input column does not match first element of input row." <<
		"\n\nColumn wins diagonal conflict.\n";
	}

    Matrix A(m,n);

    for(size_t i = 0; i< m; ++i){
        for(size_t j = 0; j< n; ++j){
            if(i==j){
                A[i][j]=c[0];
            }else if(i < j){
                A[i][j] = r[j - i ];
            }else{ // i >j
                A[i][j] = c[i - j];
            }
        }
    }    
    return A;
}
}}}