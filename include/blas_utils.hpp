// blas_utils.hpp :  include file containing templated utilities to work with vectors and matrices
//
// Copyright (C) 2017-2018 Stillwater Supercomputing, Inc.
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.

namespace sw {
	namespace hprblas {

		// These functions print matrices and vectors in a nice format
		template<typename Ty>
		void printMatrix(std::ostream& ostr, const std::string& name, const std::vector<Ty>& m) {
			size_t d = size_t(std::sqrt(m.size()));
			ostr << "Matrix: " << name << " is " << d << "x" << d << std::endl;
			ostr << std::setprecision(17);
			for (size_t i = 0; i<d; ++i) {
				for (size_t j = 0; j<d; ++j) std::cout << std::setw(20) << m[i*d + j] << " ";
				ostr << std::endl;
			}
			ostr << std::setprecision(5);
		}

		template<typename Ty>
		void printVector(std::ostream& ostr, const std::string& name, const std::vector<Ty>& v) {
			size_t d = v.size();
			ostr << "Vector: " << name << " is of size " << d << " elements" << std::endl;
			ostr << std::setprecision(17);
			for (size_t j = 0; j<d; ++j) std::cout << std::setw(20) << v[j] << " ";
			ostr << std::setprecision(5) << std::endl;
		}

	} // namespace blas

} // namespace sw
