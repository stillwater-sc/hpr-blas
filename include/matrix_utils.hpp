#pragma once
// matrix_utils.hpp :  include file containing templated utilities to work with matrices
//
// Copyright (C) 2017-2018 Stillwater Supercomputing, Inc.
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.

#include <random>
#include <algorithm>
#include "blas_utils.hpp"

namespace sw {
	namespace hprblas {

		// fill a dense MTL matrix with random values between [lowerbound, upperbound]
		template <typename Matrix>
		void uniform_rand_sorted(Matrix& A, double lowerbound = 0.0, double upperbound = 1.0)
		{
			// Use random_device to generate a seed for Mersenne twister engine.
			std::random_device rd{};
			// Use Mersenne twister engine to generate pseudo-random numbers.
			std::mt19937 engine{ rd() };
			// "Filter" MT engine's output to generate pseudo-random double values,
			// **uniformly distributed** on the closed interval [lowerbound, upperbound].
			// (Note that the range is [inclusive, inclusive].)
			std::uniform_real_distribution<double> dist{ lowerbound, upperbound };
			// Pattern to generate pseudo-random number.
			// double rnd_value = dist(engine);

			typedef typename mtl::Collection<Matrix>::value_type    value_type;
			typedef typename mtl::Collection<Matrix>::size_type     size_type;

			// generate a good set of randoms
			std::vector<value_type> v(size(A));
			for (size_type r = 0; r < num_rows(A); ++r) {
				for (size_type c = 0; c < num_cols(A); ++c) {
					v.push_back(value_type(dist(engine)));
				}
			}
			// sort them so that we have better control over the scale of each element in a row vector
			sort(v.begin(), v.end(), std::greater<value_type>());

			// for each row minus the last column, calculate the sum of elements without rounding
			sw::unum::posit<value_type::nbits, value_type::es> one(1), p;
			for (size_type r = 0; r < num_rows(A); ++r) {
				sw::unum::quire<value_type::nbits, value_type::es> q1, q2;
				size_type lastElement = num_cols(A) - 1;
				for (size_type c = 0; c < lastElement; ++c) {
					q1 += sw::unum::quire_mul(one, v[r*num_cols(A) + c]);
				}
				// truncate the value in the quire
				convert(q1.to_value(), p);
				// calculate the difference between the truncated and the non-truncated quire values
				q2 = p;
//				std::cout << "q1 :" << q1 << std::endl;
//				std::cout << "q2 :" << q2 << std::endl;
				q2 -= q1;
				convert(q2.to_value(), p);
//				std::cout << "Residual is: " << double(p) << std::endl;
				if (p.iszero()) p = one;
				v[r*num_cols(A) + lastElement] = p;
			}


			// inserters add to the elements, so we need to set the value to 0 before we begin
			A = 0.0;
			{ // extra block unfortunately needed for VS2013
			  // Create inserter for matrix m
				mtl::mat::inserter<Matrix> ins(A, num_cols(A));

				// insert sorted values in A
				size_type i = 0;
				for (size_type r = 0; r < num_rows(A); r++) {
					for (size_type c = 0; c < num_cols(A); c++) {
						ins[r][c] << v[i++];
					}
				}
				// Destructor of ins sets final state of m
			}
		}


		template<typename Matrix>
		void uniform_rand_diagonally_dominant(Matrix& A, double lowerbound = 0.0, double upperbound = 1.0) {
			// generate off-diagonal entries, calculate the sum, scale to upperbound, and set diagonal entry to slightly larger

			// Use random_device to generate a seed for Mersenne twister engine.
			std::random_device rd{};
			// Use Mersenne twister engine to generate pseudo-random numbers.
			std::mt19937 engine{ rd() };
			// "Filter" MT engine's output to generate pseudo-random double values,
			// **uniformly distributed** on the closed interval [lowerbound, upperbound].
			// (Note that the range is [inclusive, inclusive].)
			std::uniform_real_distribution<double> dist{ lowerbound, upperbound };
			// Pattern to generate pseudo-random number.
			// double rnd_value = dist(engine);

			typedef typename mtl::Collection<Matrix>::value_type    value_type;
			typedef typename mtl::Collection<Matrix>::size_type     size_type;

			// no need to null A as each element is assigned explicitely
			for (size_type r = 0; r < num_rows(A); ++r) {
				// generate a random vector of N
				size_type N = num_cols(A);
				mtl::dense_vector<value_type> v(N);
				for (size_type c = 0; c < N; ++c) {
					v[c] = value_type(dist(engine));
				}
				// add the sum of the other row elements to the diagonal
				for (size_type c = 0; c < N; ++c) {
					if (r != c) v[r] += v[c];
				}
				value_type factor = v[r] / upperbound;
				for (size_type c = 0; c < N; ++c) {
					A[r][c] = v[c] / factor;
				}
			}
		}

/*
C.3: Generating a Random Matrix with Specified Eigenvalues
Generate random orthogonal matrix G. W. Stewart (1980).

		start RandOrthog(n);
		A = I(n); // identity matrix 
		d = j(n, 1, 0);
		d[n] = sgn(RndNormal(1, 1)); // +/- 1 
		do k = n - 1 to 1 by - 1;
			// generate random Householder transformation 
			x = RndNormal(n - k + 1, 1); // column vector from N(0,1) 
			s = sqrt(x[##]); // norm(x) 
			sgn = sgn(x[1]);
			s = sgn*s;
			d[k] = -sgn;
			x[1] = x[1] + s;
			beta = s*x[1];
			// apply the Householder transformation to A 
			y = x`*A[k:n, ];
			A[k:n, ] = A[k:n, ] - x*(y / beta);
		end;
		A = d # A; // change signs of i_th row when d[i]=-1 
		return(A);
		finish;

		// helper functions 
		// return matrix of same size as A with
		// m[i,j]= { 1 if A[i,j]>=0
		//         { -1 if A[i,j]< 0
		// Similar to the SIGN function, except SIGN(0)=0 
		start sgn(A);
			return(choose(A >= 0, 1, -1));
		finish;

		// return (r x c) matrix of standard normal variates 
		start RndNormal(r, c);
			x = j(r, c);
			call randgen(x, "Normal");
			return(x);
		finish;

		// The following statements call the RANDORTHOG function to generate a random 4x4 orthogonal matrix.
		call randseed(1);
		Q = RandOrthog(4);
		print(Q);
*/
		template<typename Matrix>
		void random_orthogonal(Matrix& A) {

		}

	} // namespace blas

} // namespace sw

namespace mtl {

	namespace mat {

		// fill a dense MTL matrix with random values between [lowerbound, upperbound]
		template <typename Matrix>
		void uniform_rand(Matrix& A, double lowerbound = 0.0, double upperbound = 1.0)
		{
			// Use random_device to generate a seed for Mersenne twister engine.
			std::random_device rd{};
			// Use Mersenne twister engine to generate pseudo-random numbers.
			std::mt19937 engine{ rd() };
			// "Filter" MT engine's output to generate pseudo-random double values,
			// **uniformly distributed** on the closed interval [lowerbound, upperbound].
			// (Note that the range is [inclusive, inclusive].)
			std::uniform_real_distribution<double> dist{ lowerbound, upperbound };
			// Pattern to generate pseudo-random number.
			// double rnd_value = dist(engine);

			typedef typename Collection<Matrix>::value_type    value_type;
			typedef typename Collection<Matrix>::size_type     size_type;

			// inserters add to the elements, so we need to set the value to 0 before we begin
			A = 0.0;
			{ // extra block unfortunately needed for VS2013
				// Create inserter for matrix m
				inserter<Matrix> ins(A, num_cols(A));

				// generate and insert random values in A
				for (size_type r = 0; r < num_rows(A); r++) {
					for (size_type c = 0; c < num_cols(A); c++) {
						ins[r][c] << value_type(dist(engine));
					}
				}
				// Destructor of ins sets final state of m
			}
		}
	}
}
