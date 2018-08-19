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

			std::vector<value_type> v(size(A));
			for (size_type r = 0; r < num_rows(A); ++r) {
				for (size_type c = 0; c < num_cols(A); ++c) {
					v.push_back(value_type(dist(engine)));
				}
			}
			sort(v.begin(), v.end(), greater<>());

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
		void project_onto_number_system_accuracy(Matrix& A) {
			
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