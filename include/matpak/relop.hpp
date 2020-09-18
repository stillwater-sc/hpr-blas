#pragma once
// relop.hpp : A > B (determines if A > B element-wise)
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

namespace sw { namespace hprblas { namespace matpak {

template<typename Matrix>
bool validate_shape(const Matrix&A, const Matrix&B){
    typedef typename Matrix::size_type  size_type;
    size_type ra = num_rows(A);
	size_type ca = num_cols(A);

	size_type rb = num_rows(B);
	size_type cb = num_cols(B);

	if (ra!=rb || ca!=cb){
        return  false;
    }
    return true;
}
enum RelOp{RELOP_EQ, RELOP_NEQ, RELOP_GT, RELOP_GTE, RELOP_LT, RELOP_LTE};


template<typename Matrix>
Matrix relop(const Matrix&A, const RelOp& op, const Matrix&B){

	if(!validate_shape(A,B)) return Matrix {0};

    typedef typename Matrix::size_type  size_type;

    size_type m = num_rows(A);
    size_type n = num_cols(A);

    Matrix C(m,n);

	for(size_type i = 0; i < m; ++i){
	   for(size_type j = 0; j < n; ++j){
	      	switch (op){
            case RELOP_EQ:
                  if(A[i][j] == B[i][j]) C[i][j] = 1;
            case RELOP_NEQ:
                  if(A[i][j] != B[i][j]) C[i][j] = 1;
            case RELOP_LT:
                  if(A[i][j] < B[i][j]) C[i][j] = 1;
            case RELOP_LTE:
                  if(A[i][j] <= B[i][j]) C[i][j] = 1;
            case RELOP_GT:
                  if(A[i][j] > B[i][j]) C[i][j] = 1;
            case RELOP_GTE:
                  if(A[i][j] >= B[i][j]) C[i][j] = 1;           
            }              
	   }
    }
    return C;
}

}}}