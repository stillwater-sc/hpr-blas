// lu_decomposition.cpp example program comparing float vs posit LU decomposition algorithms
//
// Copyright (C) 2017-2019 Stillwater Supercomputing, Inc.
//
// This file is part of the HPR-BLAS project, which is released under an MIT Open Source license.

#include <chrono>
// configure the posit number system behavior
#define POSIT_ROUNDING_ERROR_FREE_IO_FORMAT 1
#include <hprblas>
#include <mtl_extensions.hpp>
#include <matrix_utils.hpp>
#include <print_utils.hpp>

template<size_t nbits, size_t es, size_t capacity = 10>
void CroutCycle(mtl::dense2D< sw::unum::posit<nbits, es> >& A, mtl::dense_vector< sw::unum::posit<nbits, es> >& x, mtl::dense_vector< sw::unum::posit<nbits, es> >& b)
{
	using namespace sw::hprblas;

	assert(num_cols(A) == size(x));
	size_t N = size(x);
	mtl::dense2D< sw::unum::posit<nbits, es> > LU(N, N);

	std::cout << "----------------- Crout cycle ------------------------\n";
	using namespace std::chrono;
	steady_clock::time_point t1 = steady_clock::now();
	Crout(A, LU);
	steady_clock::time_point t2 = steady_clock::now();
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	double elapsed = time_span.count();
	std::cout << "Crout took " << elapsed << " seconds." << std::endl;
	std::cout << "Performance " << (uint32_t)(N*N*N / (1000 * elapsed)) << " KOPS/s" << std::endl;
	SolveCrout(LU, b, x);
	printMatrix(std::cout, "Crout LU", LU);
	printVector(std::cout, "Crout Solution", x);
}

template<size_t nbits, size_t es, size_t capacity = 10>
void CroutFDPCycle(mtl::dense2D< sw::unum::posit<nbits, es> >& A, mtl::dense_vector< sw::unum::posit<nbits, es> >& x, mtl::dense_vector< sw::unum::posit<nbits, es> >& b)
{
	using namespace sw::hprblas;

	size_t d = size(b);
	assert(size(A) == d*d);
	mtl::dense2D< sw::unum::posit<nbits, es> > LU(d, d);

	std::cout << "----------------- Crout FDP cycle --------------------\n";
	using namespace std::chrono;
	steady_clock::time_point t1 = steady_clock::now();
	CroutFDP(A, LU);
	steady_clock::time_point t2 = steady_clock::now();
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	double elapsed = time_span.count();
	std::cout << "Crout with FDP took " << elapsed << " seconds." << std::endl;
	std::cout << "Performance " << (uint32_t)(d*d*d / (1000 * elapsed)) << " KOPS/s" << std::endl;
	SolveCroutFDP(LU, b, x);
	printMatrix(std::cout, "Crout FDP LU", LU);
	printVector(std::cout, "Crout FDP Solution", x);
}

template<size_t nbits, size_t es, size_t capacity = 10>
void ComparePositDecompositions(std::vector< sw::unum::posit<nbits, es> >& A, std::vector< sw::unum::posit<nbits, es> >& x, std::vector< sw::unum::posit<nbits, es> >& b) {
	size_t d = b.size();
	assert(A.size() == d*d);
	using namespace sw::hprblas;
	std::vector< sw::unum::posit<nbits, es> > LU(d*d);

	{
		using namespace std::chrono;
		steady_clock::time_point t1 = steady_clock::now();
		Crout(A, LU);
		steady_clock::time_point t2 = steady_clock::now();
		duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
		double elapsed = time_span.count();
		std::cout << "Crout took " << elapsed << " seconds." << std::endl;
		std::cout << "Performance " << (uint32_t)(d*d*d / (1000 * elapsed)) << " KOPS/s" << std::endl;

		SolveCrout(LU, b, x);
		printMatrix(std::cout, "Crout LU", LU);
		printVector(std::cout, "Solution", x);
	}

	std::cout << std::endl;
#if 0
	{
		using namespace std::chrono;
		steady_clock::time_point t1 = steady_clock::now();
		DoolittleFDP(A, LU);
		steady_clock::time_point t2 = steady_clock::now();
		duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
		double elapsed = time_span.count();
		std::cout << "Doolittle took " << elapsed << " seconds." << std::endl;
		std::cout << "Performance " << (uint32_t)(d*d*d / (1000 * elapsed)) << " KOPS/s" << std::endl;
		SolveDoolittle(LU, b, x);
		printMatrix(std::cout, "Doolittle LU", LU);
		printVector(std::cout, "Solution", x);
	}


	std::cout << std::endl;

	{
		using namespace std::chrono;
		steady_clock::time_point t1 = steady_clock::now();
		CholeskyFDP(A, LU);
		steady_clock::time_point t2 = steady_clock::now();
		duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
		double elapsed = time_span.count();
		std::cout << "Cholesky took " << elapsed << " seconds." << std::endl;
		std::cout << "Performance " << (uint32_t)(d*d*d / (1000 * elapsed)) << " KOPS/s" << std::endl;
		SolveCholesky(LU, b, x);
		printMatrix(std::cout, "Cholesky LU", LU);
		printVector(std::cout, "Solution", x);
	}
#endif
}

template<typename Scalar>
void RandomMatrix() {
	mtl::dense_vector<Scalar> x(5), b(5), xprime(5);
	mtl::dense2D<Scalar> A(5, 5);
	mtl::mat::uniform_rand(A, -1.0, 1.0);

	x = 1.0;
	b = A * x;
	cout << endl;
	printMatrix(cout, "Matrix A(5x5):\n", A);
	cout << endl;
	cout << endl;
	printVector(cout, "RHS    b(5)  :\n", b);
	cout << endl;
	CroutCycle<nbits, es, capacity>(A, xprime, b);
	printVector(cout, "RHS    x(5)  :\n", xprime);
	cout << endl;
	CroutFDPCycle<nbits, es, capacity>(A, xprime, b);
	printVector(cout, "RHS    x(5)  :\n", xprime);
}

template<typename Ty>
void CompareIEEEDecompositions(std::vector<Ty>& A, std::vector<Ty>& x, std::vector<Ty>& b) {
	size_t d = b.size();
	assert(A.size() == d*d);
	using namespace sw::hprblas;
	std::vector<Ty> LU(d*d);

	{
		using namespace std::chrono;
		steady_clock::time_point t1 = steady_clock::now();
		Crout(A, LU);
		steady_clock::time_point t2 = steady_clock::now();
		duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
		double elapsed = time_span.count();
		std::cout << "Crout took " << elapsed << " seconds." << std::endl;
		std::cout << "Performance " << (uint32_t)(d*d*d / (1000 * elapsed)) << " KOPS/s" << std::endl;

		SolveCrout(LU, b, x);
		printMatrix(std::cout, "Crout LU", LU);
		printVector(std::cout, "Solution", x);
	}


	std::cout << std::endl;

	{
		using namespace std::chrono;
		steady_clock::time_point t1 = steady_clock::now();
		Doolittle(A, LU);
		steady_clock::time_point t2 = steady_clock::now();
		duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
		double elapsed = time_span.count();
		std::cout << "Doolittle took " << elapsed << " seconds." << std::endl;
		std::cout << "Performance " << (uint32_t)(d*d*d / (1000 * elapsed)) << " KOPS/s" << std::endl;
		SolveCrout(LU, b, x);
		printMatrix(std::cout, "Doolittle LU", LU);
		printVector(std::cout, "Solution", x);

		SolveDoolittle(LU, b, x);
		printMatrix(std::cout, "Doolittle LU", LU);
		printVector(std::cout, "Solution", x);
	}


	std::cout << std::endl;

	{
		using namespace std::chrono;
		steady_clock::time_point t1 = steady_clock::now();
		Cholesky(A, LU);
		steady_clock::time_point t2 = steady_clock::now();
		duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
		double elapsed = time_span.count();
		std::cout << "Cholesky took " << elapsed << " seconds." << std::endl;
		std::cout << "Performance " << (uint32_t)(d*d*d / (1000 * elapsed)) << " KOPS/s" << std::endl;
		SolveCrout(LU, b, x);
		printMatrix(std::cout, "Cholesky LU", LU);
		printVector(std::cout, "Solution", x);

		SolveCholesky(LU, b, x);
		printMatrix(std::cout, "Cholesky LU", LU);
		printVector(std::cout, "Solution", x);
	}
}


int main(int argc, char** argv)
try {
	using namespace std;
	using namespace mtl;
	using namespace sw::unum;
	using namespace sw::hprblas;

	// a 32-bit float and a <27,1> posit have the same number of significand bits around 1.0
	constexpr size_t nbits = 16;
	constexpr size_t es = 1;
	constexpr size_t capacity = 10;

	{
		using Scalar = posit<nbits, es>;
		size_t N = 5;
		dense2D<Scalar> U(N, N), L(N, N), A(N, N);
		fill_U(U);
		fill_L(L);

		// We want to solve the system Ax=b
		matmul(A, L, U);   // construct the A matrix to solve
		printMatrix(cout, "A = LU", A);
		cout << endl;

		// define a difficult solution
		Scalar eps = std::numeric_limits<Scalar>::epsilon();
		Scalar epsminus = Scalar(1.0) - eps;
		Scalar epsplus = Scalar(1.0) + eps;
		dense_vector<Scalar> x(N), b(N);
		x = epsplus;
		matvec(A, x, b);   // construct the right hand side
		printVector(cout, "x", x);
		printVector(cout, "b", b);
		cout << endl;
		dense_vector<Scalar> xprime(N);
		CroutCycle(A, xprime, b);
		printVector(cout, "xprime", xprime);
		cout << endl;
		CroutFDPCycle(A, xprime, b);
		printVector(cout, "xprime", xprime);
	}

#if 0
	cout << "LinearSolve regular dot product" << endl;
	CompareIEEEDecompositions(Aieee, xieee, bieee);
	cout << endl << ">>>>>>>>>>>>>>>>" << endl;
	cout << "LinearSolve fused-dot product" << endl;
	ComparePositDecompositions(Aposit, xposit, bposit);
#endif


	return EXIT_SUCCESS;
}
catch (char const* msg) {
	std::cerr << msg << std::endl;
	return EXIT_FAILURE;
}
catch (const posit_arithmetic_exception& err) {
	std::cerr << "Uncaught posit arithmetic exception: " << err.what() << std::endl;
	return EXIT_FAILURE;
}
catch (const quire_exception& err) {
	std::cerr << "Uncaught quire exception: " << err.what() << std::endl;
	return EXIT_FAILURE;
}
catch (const posit_internal_exception& err) {
	std::cerr << "Uncaught posit internal exception: " << err.what() << std::endl;
	return EXIT_FAILURE;
}
catch (std::runtime_error& err) {
	std::cerr << err.what() << std::endl;
	return EXIT_FAILURE;
}
catch (...) {
	std::cerr << "Caught unknown exception" << std::endl;
	return EXIT_FAILURE;
}
