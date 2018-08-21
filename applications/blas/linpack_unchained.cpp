// linpack_unchained.cpp: system solver benchmark with exact solution
// This program completes the design goal of LINPACK, that is,
// create an idealized problem where both the input values and the correct answer 
// are expressible in the numerical vocabulary of the computing environment.
//
// Copyright (C) 2017-18 Stillwater Supercomputing, Inc.
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.

#include "common.hpp"

#include <iostream>
#include <typeinfo>

#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>
#include <posit>
#include "blas.hpp"
#include "vector_utils.hpp"
#include "matrix_utils.hpp"

using namespace std;
using namespace sw::unum;

#define BACKEND_MTL 0
#define BACKEND_EIGEN 1

#if defined(BACKEND_MTL)
#if ARITHMETIC_POSIT
using Tensor = mtl::Tensor<  posit<nbits, es> >;
#elif ARITHMETIC_INT8
using Tensor = mtl::Tensor< uint8_t >;
#elif ARITHMETIC_FP16
using Tensor = mtl::Tensor< fp16 >;
#endif
#elif defined(BACKEND_EIGEN)
#if ARITHMETIC_POSIT
using Tensor = Eigen::Tensor<  posit<nbits, es> >;
#elif ARITHMETIC_INT8
using Tensor = Eigen::Tensor< uint8_t >;
#elif ARITHMETIC_FP16
using Tensor = Eigen::Tensor< fp16 >;
#endif
#endif

int main(int argc, char** argv)
try {
	const size_t nbits = 16;
	const size_t es = 1;
	const size_t vecSize = 32;

	{
		using Scalar = posit<nbits, es>;
		using Vector = mtl::dense_vector< Scalar >;
		using Matrix = mtl::dense2D< Scalar >;

		constexpr size_t N = 10;

		// generate a uniform random matrix
		Matrix  A(N, N);
		sw::hprblas::uniform_rand_sorted(A);
		Vector r(N), rtilde(N);
		Vector	x(N, 1.0), b(N);
		b = A * x;
		x = 0.99;
		r = b - A * x;
		rtilde = r;
		Scalar rho = sw::hprblas::fused_dot<Vector, nbits, es>(rtilde, r);
		cout << "rho: " << double(rho) << endl;
		sw::hprblas::printVector(cout, "r: ", r);

		return 0;

#if ITL_FIXED
		// Create an ILU(0) preconditioner
		itl::pc::ilu_0<Matrix>        P(A);

		// Set b such that x == 1 is solution; start with x == 0
		Vector	x(N, 1.0), b(N);
		b = A * x; x = 0;

		// Termination criterion: r < 1e-6 * b or N iterations
		itl::noisy_iteration<double>       iter(b, 500, 1.e-6);

		// Solve Ax == b with left preconditioner P
		itl::bicgstab(A, x, b, P, iter);
#endif

		if (b != x) {
			typedef typename mtl::Collection<Matrix>::size_type     size_type;
			posit<nbits, es> p, one(1);
			for (size_type r = 0; r < num_rows(A); ++r) {
				sw::unum::quire<nbits, es> q, qt;
				for (size_type c = 0; c < num_cols(A); ++c) {
					p = A[r][c];
					q += quire_mul(one, p);
					qt.reset();
					qt += quire_mul(one, p);
					cout << qt << endl;
				}
				cout << q << endl;

				convert(q.to_value(), p);
				qt.reset();
				qt += quire_mul(posit<nbits, es>(1.0), p);
				cout << qt << endl << endl;
			}
		}
		else {
			sw::hprblas::printVector(cout, "Solution vector: ", x);
		}
	}

	return 0;

	{
		using Matrix = mtl::dense2D<float>;
		using Vector = mtl::dense_vector<float>;
		Matrix  A(4, 4), L(4, 4), U(4, 4), LU(4, 4);
		Vector	x(4), b(4), xx(4);
		double 	c = 1.0;

		mtl::mat::uniform_rand(A);  // uniform random with values between [0,1]
		LU = A;
		lu(LU);
		cout << A << endl;
		cout << LU << endl;
	}


	int nrOfFailedTestCases = 0;
	return nrOfFailedTestCases > 0 ? EXIT_FAILURE : EXIT_SUCCESS;
}
catch (char const* msg) {
	cerr << msg << endl;
	return EXIT_FAILURE;
}
catch (const std::runtime_error& err) {
	cerr << err.what() << endl;
	return EXIT_FAILURE;
}
catch (...) {
	cerr << "caught unknown exception" << endl;
	return EXIT_FAILURE;
}
