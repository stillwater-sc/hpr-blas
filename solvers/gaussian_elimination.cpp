// gaussian_elimination.cpp: example program comparing float vs posit Gaussian Elimination (LU Decomposition) equation solver
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

#include <chrono>
#ifdef _WINDOWS
#define MTL_WITH_INITLIST
#endif
// enable posit arithmetic exceptions
#define POSIT_THROW_ARITHMETIC_EXCEPTION 1
#include <hprblas>
// matrix generators
#include <generators/matrix_generators.hpp>
// utilities
#include <universal/functions/isrepresentable.hpp>
#include <utils/print_utils.hpp>

template<size_t nbits, size_t es, size_t capacity = 10>
void ComparePositDecompositions(mtl::mat::dense2D< sw::unum::posit<nbits, es> >& A, mtl::vec::dense_vector< sw::unum::posit<nbits, es> >& x, mtl::vec::dense_vector< sw::unum::posit<nbits, es> >& b) {
	assert(mtl::mat::num_rows(A) == mtl::mat::num_cols(A));

	using namespace sw::hprblas;
	size_t N = mtl::mat::num_cols(A);
	mtl::mat::dense2D< sw::unum::posit<nbits, es> > LU(N,N);

	{
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
		CholeskyFDP(A, LU);
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
#endif
}


template<typename Matrix, typename Vector>
void CompareIEEEDecompositions(Matrix& A, Vector& x, Vector& b) {
	using namespace sw::hprblas;
	assert(mtl::mat::num_rows(A) == mtl::mat::num_cols(A));
	size_t N = A.num_cols();
	Matrix LU(N,N);

	{
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
		printVector(std::cout, "Solution", x);
	}


	std::cout << std::endl;
#if 0
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


#endif
}

/* TBD
template<> class std::numeric_limits< sw::unum::posit<28, 1> >
	: public _Num_float_base
{	// limits for type posit<28, 1>
public:
	typedef sw::unum::posit<28,1> _Ty;

	static constexpr _Ty(min)() _THROW0()
	{	// return minimum value
		return (sw::unum::minpos<28,1>());
	}

	static constexpr _Ty(max)() _THROW0()
	{	// return maximum value
		return (_DBL_MAX);
	}

	static constexpr _Ty lowest() _THROW0()
	{	// return most negative value
		return (-(max)());
	}

	static constexpr _Ty epsilon() _THROW0()
	{	// return smallest effective increment from 1.0
		return (_DBL_EPSILON);
	}

	static constexpr _Ty round_error() _THROW0()
	{	// return largest rounding error
		return (0.5);
	}

	static constexpr _Ty denorm_min() _THROW0()
	{	// return minimum denormalized value
		return (_DBL_TRUE_MIN);
	}

	static constexpr _Ty infinity() _THROW0()
	{	// return positive infinity
		return (__builtin_huge_val());
	}

	static constexpr _Ty quiet_NaN() _THROW0()
	{	// return non-signaling NaN
		return (__builtin_nan("0"));
	}

	static constexpr _Ty signaling_NaN() _THROW0()
	{	// return signaling NaN
		return (__builtin_nans("1"));
	}

	_STCONS(int, digits, DBL_MANT_DIG);
	_STCONS(int, digits10, DBL_DIG);

	_STCONS(int, max_digits10, 2 + DBL_MANT_DIG * 301L / 1000);

	_STCONS(int, max_exponent, (int)DBL_MAX_EXP);
	_STCONS(int, max_exponent10, (int)DBL_MAX_10_EXP);
	_STCONS(int, min_exponent, (int)DBL_MIN_EXP);
	_STCONS(int, min_exponent10, (int)DBL_MIN_10_EXP);
};
*/

int main(int argc, char** argv)
try {
	using namespace std;
	using namespace sw::unum;
	using namespace sw::hprblas;

	using IEEEType = float;
	float eps = std::numeric_limits<float>::epsilon();
	float epsminus = 1.0f - eps;
	float epsplus  = 1.0f + eps;

	// We want to solve the system Ax=b
	int N = 5;

	{
		// IEEE input data set up
		cout << "IEEE inputs\n";
		mtl::mat::dense2D<IEEEType> Uieee = {     // define the upper triangular matrix
			{ 1.0, 2.0, 3.0, 4.0, 5.0 },
			{ 0.0, 1.0, 2.0, 3.0, 4.0 },
			{ 0.0, 0.0, 1.0, 2.0, 3.0 },
			{ 0.0, 0.0, 0.0, 1.0, 2.0 },
			{ 0.0, 0.0, 0.0, 0.0, 1.0 },
		};
		mtl::mat::dense2D<IEEEType> Lieee = {     // define the lower triangular matrix
			{ 1.0, 0.0, 0.0, 0.0, 0.0 },
			{ 2.0, 1.0, 0.0, 0.0, 0.0 },
			{ 3.0, 2.0, 1.0, 0.0, 0.0 },
			{ 4.0, 3.0, 2.0, 1.0, 0.0 },
			{ 5.0, 4.0, 3.0, 2.0, 1.0 },
		};
		mtl::mat::dense2D<IEEEType> Aieee(N,N);
		Aieee = Lieee * Uieee;   // construct the A matrix to solve
		printMatrix(cout, "L", Lieee);
		printMatrix(cout, "U", Uieee);
		printMatrix(cout, "A", Aieee);

		// define a difficult solution
		mtl::vec::dense_vector<IEEEType> xieee = {
			epsplus,
			epsplus,
			epsplus,
			epsplus,
			epsplus
		};
		mtl::vec::dense_vector<IEEEType> bieee(N);
		bieee = Aieee * xieee;   // construct the right hand side
		printVector(cout, "b", bieee);
		cout << "LinearSolve regular dot product" << endl;
		CompareIEEEDecompositions(Aieee, xieee, bieee);
	}

	{  // Posit comparison
		// a 32-bit float and a <27,1> posit have the same number of significand bits around 1.0
		constexpr size_t nbits = 32;
		constexpr size_t es = 2;
		//constexpr size_t capacity = 10;
		using PositType = posit<nbits, es>;
		cout << "Using " << dynamic_range(posit<nbits, es>()) << endl;

		// repeat set up for posits
		cout << "Posit inputs\n";
		mtl::mat::dense2D<PositType> Uposit = {     // define the upper triangular matrix
			{ 1.0, 2.0, 3.0, 4.0, 5.0 },
			{ 0.0, 1.0, 2.0, 3.0, 4.0 },
			{ 0.0, 0.0, 1.0, 2.0, 3.0 },
			{ 0.0, 0.0, 0.0, 1.0, 2.0 },
			{ 0.0, 0.0, 0.0, 0.0, 1.0 },
		};
		mtl::mat::dense2D<PositType> Lposit = {     // define the lower triangular matrix
			{ 1.0, 0.0, 0.0, 0.0, 0.0 },
			{ 2.0, 1.0, 0.0, 0.0, 0.0 },
			{ 3.0, 2.0, 1.0, 0.0, 0.0 },
			{ 4.0, 3.0, 2.0, 1.0, 0.0 },
			{ 5.0, 4.0, 3.0, 2.0, 1.0 },
		};
		mtl::mat::dense2D<PositType> Aposit(N,N);
		Aposit = fmm(Lposit, Uposit);   // construct the A matrix to solve
		printMatrix(cout, "L", Lposit);
		printMatrix(cout, "U", Uposit);
		printMatrix(cout, "A", Aposit);
		// define a difficult solution
		mtl::vec::dense_vector<PositType> xposit = {
			epsplus,
			epsplus,
			epsplus,
			epsplus,
			epsplus
		};
		mtl::vec::dense_vector<PositType> bposit(N);
		bposit = fmv(Aposit, xposit);   // construct the right hand side
		printVector(cout, "b", bposit);
		cout << endl << ">>>>>>>>>>>>>>>>" << endl;
		cout << "LinearSolve fused-dot product" << endl;
		ComparePositDecompositions(Aposit, xposit, bposit);
	}

#if 1
	cout << "posit<25,1>\n";
	cout << "1.0 - FLT_EPSILON = " << setprecision(17) << epsminus << " converts to " << posit<25, 1>(epsminus) << endl;
	cout << "1.0 + FLT_EPSILON = " << setprecision(17) << epsplus << " converts to " << posit<25, 1>(epsplus) << endl;
	cout << "posit<26,1>\n";
	cout << "1.0 - FLT_EPSILON = " << setprecision(17) << epsminus << " converts to " << posit<26, 1>(epsminus) << endl;
	cout << "1.0 + FLT_EPSILON = " << setprecision(17) << epsplus << " converts to " << posit<26, 1>(epsplus) << endl;
	cout << "posit<27,1>\n";
	cout << "1.0 - FLT_EPSILON = " << setprecision(17) << epsminus << " converts to " << posit<27, 1>(epsminus) << endl;
	cout << "1.0 + FLT_EPSILON = " << setprecision(17) << epsplus << " converts to " << posit<27, 1>(epsplus) << endl;
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
