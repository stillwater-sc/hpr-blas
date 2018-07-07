// mm.cpp example program comparing float vs posit matrix multiply algorithms
//
// Copyright (C) 2017-2018 Stillwater Supercomputing, Inc.
//
// This file is part of the HPR-BLAS project, which is released under an MIT Open Source license.
#include "common.hpp"
#include <hprblas>

// fill a dense MTL matrix with random values between (0, 1)
template <typename Ty>
void rand(mtl::dense2D<Ty>& m)
{
	// Use random_device to generate a seed for Mersenne twister engine.
	std::random_device rd{};
	// Use Mersenne twister engine to generate pseudo-random numbers.
	std::mt19937 engine{ rd() };
	// "Filter" MT engine's output to generate pseudo-random double values,
	// **uniformly distributed** on the closed interval [0, 1].
	// (Note that the range is [inclusive, inclusive].)
	std::uniform_real_distribution<double> dist{ 0.0, 1.0 };
	// Pattern to generate pseudo-random number.
	// double rnd_value = dist(engine);

	// inserters add to the elements, so we need to set the value to 0 before we begin
	m = 0.0f;
	// Create inserter for matrix m
	mtl::mat::inserter< mtl::dense2D<Ty> > ins(m);

	// generate and insert random values in m
	for (int i = 0; i < m.num_rows(); ++i) {
		for (int j = 0; j < m.num_cols(); ++j) {
			ins[i][j] << float(dist(engine));
		}
	}

	// Destructor of ins sets final state of m
}

template <typename Matrix>
void fill_and_print(Matrix& A, char name)
{
	// Set values in traditional way
	A = 1.2, 3.4,
		5.6, 7.8;

	// Just print them
	std::cout << name << " is \n" << A << "\n";
}

int matrix_types(int, char**)
{
#if defined(MTL_WITH_VARIADIC_TEMPLATE) && defined(MTL_WITH_TEMPLATE_ALIAS)
	using namespace mtl;

	// Compressed matrix
	matrix<float, compressed>                 A(2, 2);
	fill_and_print(A, 'A');

	// Banded matrix
	matrix<float, sparse, banded>             B(2, 2);
	fill_and_print(B, 'B');

	// Matrix in the ELLPACK format 
	matrix<double, ellpack>                   C(2, 2);
	fill_and_print(C, 'C');

	// Coordinate matrix
	matrix<float, coordinate>                 D(2, 2);
	fill_and_print(D, 'D');

	// Morton-order matrix with the default mask
	matrix<double, morton>                    E(2, 2);
	fill_and_print(E, 'E');

	// Matrix with a Morton mask is of course a Morton-order matrix
	matrix<double, mask<shark_z_64_row_mask>> F(2, 2);
	fill_and_print(F, 'F');
#endif
	return 0;
}

int morton_recursion()
{
	using namespace mtl; 
	using namespace mtl::mat;

	const unsigned n = 20;
	dense2D<double>                            A(n, n), B(n, n);
	morton_dense<double, doppled_64_row_mask>  C(n, n);

	hessian_setup(A, 3.0); hessian_setup(B, 1.0);
	hessian_setup(C, 2.0);

	// Corresponds to A= B * B;
	mult(B, B, A);

	A = B * B;   // use BLAS
	A = B * C;   // use recursion + tiling from MTL4

	A += B * C;  // Increment A by the product of B and C
	A -= B * C;  // Likewise with decrement

	return 0;
}

int main(int argc, char** argv)
try {
	using namespace std;
	using namespace sw::unum;
	using namespace sw::hprblas;
	using namespace mtl;

	// a 32-bit float and a <27,1> posit have the same number of significand bits around 1.0
	constexpr size_t nbits = 27;
	constexpr size_t es = 2;
	constexpr size_t capacity = 10;

	typedef float            IEEEType;
	typedef posit<nbits, es> PositType;
	cout << spec_to_string(posit<nbits, es>()) << endl;

	float eps = std::numeric_limits<IEEEType>::epsilon();
	float epsminus = 1.0f - eps;
	float epsplus = 1.0f + eps;
	float max_fexp = std::numeric_limits<IEEEType>::max_exponent;
	float min_fexp = std::numeric_limits<IEEEType>::min_exponent;
	cout << "IEEE float: epsilon " << eps << " min exp " << min_fexp << " max exp " << max_fexp << endl;

	constexpr int dim = 1000;
	using namespace std::chrono;
	steady_clock::time_point t1 = steady_clock::now();
	mtl::dense2D<IEEEType> A(dim, dim), B(dim, dim), C(dim, dim);
	steady_clock::time_point t2 = steady_clock::now();
	duration<float> time_span = duration_cast< duration<float> >(t2 - t1);
	float elapsed = time_span.count();
	cout << "Construction took " << elapsed << " seconds.\n";

	t1 = steady_clock::now();
	rand(A);
	rand(B);
	t2 = steady_clock::now();
	time_span = duration_cast< duration<float> >(t2 - t1);
	elapsed = time_span.count();
	cout << "Random fill took " << elapsed << " seconds.\n";

	t1 = steady_clock::now();
	C = A * B;
	t2 = steady_clock::now();
	time_span = duration_cast<duration<float>> (t2 - t1);
	elapsed = time_span.count();
	float flops = dim*dim*dim / elapsed;
	cout << "Performance: " << flops << "SP FLOPS\n";

	return EXIT_SUCCESS;
}
catch (char const* msg) {
	std::cerr << msg << std::endl;
	return EXIT_FAILURE;
}
catch (...) {
	std::cerr << "Caught unknown exception" << std::endl;
	return EXIT_FAILURE;
}
