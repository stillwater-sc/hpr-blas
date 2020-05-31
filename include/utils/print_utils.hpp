#pragma once
// print_utils.hpp :  include file containing templated utilities to work with vectors and matrices
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.


// These functions print matrices and vectors in a nice format
namespace sw {
namespace hprblas {

// printVector pretty prints a vector
template<typename Vector>
void printVector(std::ostream& ostr, const std::string& name, const Vector& v) {
	size_t d = size(v);
	ostr << "Vector: " << name << " is of size " << d << " elements" << std::endl;
	std::streamsize old_prec = ostr.precision();
	ostr << std::setprecision(17);
	for (auto&& element: v) std::cout << std::setw(20) << element << " ";
	ostr << std::setprecision(old_prec) << std::endl;
}

// printMatrix pretty prints a 2D dense matrix
template<typename Matrix>
void printMatrix(std::ostream& ostr, const std::string& name, const Matrix& M, int precision = 17, bool hex = false) {
	using namespace mtl;
	size_t d = num_rows(M);
	ostr << "Matrix: " << name << " is " << d << "x" << d << std::endl;
	std::streamsize old_prec = ostr.precision();
	ostr << std::setprecision(precision);
	for (size_t i = 0; i<d; ++i) {
		if (hex) {
			std::cout << '[';
			for (size_t j = 0; j < d; ++j) std::cout << " " << sw::unum::hex_format(M[i][j]);
			std::cout << ']';
		}
		else {
			for (size_t j = 0; j<d; ++j) std::cout << std::setw(20) << M[i][j] << " ";
		}
		
		ostr << std::endl;
	}
	ostr << std::setprecision(old_prec);
}

// printMatrix pretty prints a 2D dense matrix that is represented by an STL vector
template<typename Ty>
void printMatrix(std::ostream& ostr, const std::string& name, const std::vector<Ty>& M) {
	using namespace mtl;
	size_t d = size(M);
	ostr << "Matrix: " << name << " is " << d << "x" << d << std::endl;
	std::streamsize old_prec = ostr.precision();
	ostr << std::setprecision(17);
	for (size_t i = 0; i<d; ++i) {
		for (size_t j = 0; j<d; ++j) std::cout << std::setw(20) << M[i*d + j] << " ";
		ostr << std::endl;
	}
	ostr << std::setprecision(old_prec);
}


// printSubMatrix pretty prints a sub-matrix of a 2D dense matrix
template<typename Matrix>
void printSubMatrix(std::ostream& ostr, const std::string& name, const Matrix& A, unsigned ai, unsigned aj, unsigned blockHeight, unsigned blockWidth) {
	ostr << "Sub-Matrix (" << ai << "," << aj << ") of " << name << " of size " << blockHeight << "x" << blockWidth << std::endl;
	std::streamsize old_prec = ostr.precision();
	ostr << std::setprecision(17);

	unsigned aRows = unsigned(mtl::mat::num_rows(A));
	unsigned aCols = unsigned(mtl::mat::num_cols(A));

	unsigned aRow = ai*blockHeight;
	unsigned aCol = aj*blockWidth;

	unsigned maxRow = (aRow + blockHeight < aRows) ? blockHeight : aRows - aRow;
	unsigned maxCol = (aCol + blockWidth < aCols) ? blockWidth : aCols - aCol;

	for (unsigned i = 0; i < maxRow; ++i) {
		for (unsigned j = 0; j < maxCol; ++j) {
			std::cout << std::setw(20) << A[aRow + i][aCol + j] << " ";
		}
		ostr << std::endl;
	}
	ostr << std::setprecision(old_prec);
}
	
} // namespace blas
} // namespace sw
