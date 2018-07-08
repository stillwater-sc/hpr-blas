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

	constexpr int dim = 50;
	using namespace std::chrono;
	for (int dim = 16; dim < 2049; dim *= 2) {
		cout << "Matrix dimensions are: " << dim << " x " << dim << endl;

		steady_clock::time_point t1 = steady_clock::now();
		mtl::dense2D<IEEEType> A(dim, dim), B(dim, dim), C(dim, dim);
		steady_clock::time_point t2 = steady_clock::now();
		duration<float> time_span = duration_cast< duration<float> >(t2 - t1);
		float elapsed = time_span.count();
		cout << "  Construction     " << elapsed << " seconds.\n";

		t1 = steady_clock::now();
		rand(A);
		rand(B);
		t2 = steady_clock::now();
		time_span = duration_cast< duration<float> >(t2 - t1);
		elapsed = time_span.count();
		cout << "  Random fill      " << elapsed << " seconds.\n";

		t1 = steady_clock::now();
		int N = 10;
		for (int i = 0; i < N; ++i) {
			C = A * B;
		}
		t2 = steady_clock::now();
		time_span = duration_cast<duration<float>> (t2 - t1);
		elapsed = time_span.count();
		cout << "  Matmul iteration " << elapsed << " seconds.\n";
		float flops = float(N)*dim*dim*dim / elapsed * 1.0e-9;
		cout << "  Performance:     " << flops << " GFLOPS Single Precision\n";
	}

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
