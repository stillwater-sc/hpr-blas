#pragma once
// blas_utils.hpp :  include file containing templated utilities to work with vectors and matrices
//
// Copyright (C) 2017-2018 Stillwater Supercomputing, Inc.
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.

namespace sw {
	namespace hprblas {

		// These functions print matrices and vectors in a nice format
		template<typename Matrix>
		void printMatrix(std::ostream& ostr, const std::string& name, const Matrix& M) {
			size_t d = num_rows(M);
			ostr << "Matrix: " << name << " is " << d << "x" << d << std::endl;
			std::streamsize old_prec = ostr.precision();
			ostr << std::setprecision(17);
			for (size_t i = 0; i<d; ++i) {
				for (size_t j = 0; j<d; ++j) std::cout << std::setw(20) << double(M[i][j]) << " ";
				ostr << std::endl;
			}
			ostr << std::setprecision(old_prec);
		}

		// print Matrix that is an STL vector
		template<typename Ty>
		void printMatrix(std::ostream& ostr, const std::string& name, const std::vector<Ty>& M) {
			size_t d = size(M);
			ostr << "Matrix: " << name << " is " << d << "x" << d << std::endl;
			std::streamsize old_prec = ostr.precision();
			ostr << std::setprecision(17);
			for (size_t i = 0; i<d; ++i) {
				for (size_t j = 0; j<d; ++j) std::cout << std::setw(20) << double(M[i*d + j]) << " ";
				ostr << std::endl;
			}
			ostr << std::setprecision(old_prec);
		}

		template<typename Vector>
		void printVector(std::ostream& ostr, const std::string& name, const Vector& v) {
			size_t d = size(v);
			ostr << "Vector: " << name << " is of size " << d << " elements" << std::endl;
			std::streamsize old_prec = ostr.precision();
			ostr << std::setprecision(17);
			for (size_t j = 0; j<d; ++j) std::cout << std::setw(20) << double(v[j]) << " ";
			ostr << std::setprecision(old_prec) << std::endl;
		}

	} // namespace blas

} // namespace sw
