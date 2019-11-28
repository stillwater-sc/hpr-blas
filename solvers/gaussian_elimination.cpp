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
#include <utils/print_utils.hpp>


#ifdef LOCAL
// The following compact LU factorization schemes are described
// in Dahlquist, Bjorck, Anderson 1974 "Numerical Methods".
//
// S and D are d by d matrices.  However, they are stored in
// memory as 1D arrays of length d*d.  Array indices are in
// the C style such that the first element of array A is A[0]
// and the last element is A[d*d-1].
//
// These routines are written with separate source S and
// destination D matrices so the source matrix can be retained
// if desired.  However, the compact schemes were designed to
// perform in-place computations to save memory.  In
// other words, S and D can be the SAME matrix.  

// Crout implements an in-place LU decomposition, that is, S and D can be the same
// Crout uses unit diagonals for the upper triangle
template<typename Ty>
void Crout(std::vector<Ty>& S, std::vector<Ty>& D) {
	size_t d = size_t(std::sqrt(S.size()));
	assert(S.size() == d*d);
	assert(D.size() == d*d);
	for (size_t k = 0; k < d; ++k) {
		for (size_t i = k; i < d; ++i) {
			Ty sum = 0.;
			for (size_t p = 0; p < k; ++p) sum += D[i*d + p] * D[p*d + k];
			D[i*d + k] = S[i*d + k] - sum; // not dividing by diagonals
		}
		for (size_t j = k + 1; j < d; ++j) {
			Ty sum = 0.;
			for (size_t p = 0; p < k; ++p) sum += D[k*d + p] * D[p*d + j];
			D[k*d + j] = (S[k*d + j] - sum) / D[k*d + k];
		}
	}
}

// SolveCrout takes an LU decomposition, LU, and a right hand side vector, b, and produces a result, x.
template<typename Ty>
void SolveCrout(const std::vector<Ty>& LU, const std::vector<Ty>& b, std::vector<Ty>& x) {
	int d = (int)b.size();
	std::vector<Ty> y(d);
	for (int i = 0; i < d; ++i) {
		Ty sum = 0.0;
		for (int k = 0; k < i; ++k) sum += LU[i*d + k] * y[k];
		y[i] = (b[i] - sum) / LU[i*d + i];

	}
	for (int i = d - 1; i >= 0; --i) {
		Ty sum = 0.0;
		for (int k = i + 1; k < d; ++k) {
			//cout << "lu[] = " << LU[i*d+k] << " x[" << k << "] = " << x[k] << endl;
			sum += LU[i*d + k] * x[k];
		}
		//cout << "sum " << sum << endl;
		x[i] = (y[i] - sum); // not dividing by diagonals
	}
}

template<size_t nbits, size_t es, size_t capacity = 10>
void CroutFDP(std::vector< sw::unum::posit<nbits, es> >& S, std::vector< sw::unum::posit<nbits, es> >& D) {
	size_t d = size_t(std::sqrt(S.size()));
	assert(S.size() == d*d);
	assert(D.size() == d*d);
	using namespace sw::unum;
	for (int k = 0; k < d; ++k) {
		for (int i = k; i < d; ++i) {
			quire<nbits, es, capacity> q = 0.0;
			//for (int p = 0; p < k; ++p) q += D[i*d + p] * D[p*d + k];   if we had expression templates for the quire
			for (int p = 0; p < k; ++p) q += quire_mul(D[i*d + p], D[p*d + k]);
			posit<nbits, es> sum;
			sum.convert(q.to_value());     // one and only rounding step of the fused-dot product
			D[i*d + k] = S[i*d + k] - sum; // not dividing by diagonals
		}
		for (int j = k + 1; j < d; ++j) {
			quire<nbits, es, capacity> q = 0.0;
			//for (int p = 0; p < k; ++p) q += D[k*d + p] * D[p*d + j];   if we had expression templates for the quire
			for (int p = 0; p < k; ++p) q += quire_mul(D[k*d + p], D[p*d + j]);
			posit<nbits, es> sum;
			sum.convert(q.to_value());   // one and only rounding step of the fused-dot product
			D[k*d + j] = (S[k*d + j] - sum) / D[k*d + k];
		}
	}
}

// SolveCrout takes an LU decomposition, LU, and a right hand side vector, b, and produces a result, x.
template<size_t nbits, size_t es, size_t capacity = 10>
void SolveCroutFDP(const std::vector< sw::unum::posit<nbits, es> >& LU, const std::vector< sw::unum::posit<nbits, es> >& b, std::vector< sw::unum::posit<nbits, es> >& x) {
	using namespace sw::unum;
	int d = int(b.size());
	std::vector< posit<nbits, es> > y(d);
	for (int i = 0; i < d; ++i) {
		quire<nbits, es, capacity> q = 0.0;
		// for (int k = 0; k < i; ++k) q += LU[i*d + k] * y[k];   if we had expression templates for the quire
		for (int k = 0; k < i; ++k) q += quire_mul(LU[i*d + k], y[k]);
		posit<nbits, es> sum;
		sum.convert(q.to_value());   // one and only rounding step of the fused-dot product
		y[i] = (b[i] - sum) / LU[i*d + i];
	}
	for (int i = d - 1; i >= 0; --i) {
		quire<nbits, es, capacity> q = 0.0;
		// for (int k = i + 1; k < d; ++k) q += LU[i*d + k] * x[k];   if we had expression templates for the quire
		for (int k = i + 1; k < d; ++k) {
			//cout << "lu[] = " << LU[i*d + k] << " x[" << k << "] = " << x[k] << endl;
			q += quire_mul(LU[i*d + k], x[k]);
		}
		posit<nbits, es> sum;
		sum.convert(q.to_value());   // one and only rounding step of the fused-dot product
		//cout << "sum " << sum << endl;
		x[i] = (y[i] - sum); // not dividing by diagonals
	}
}

// Doolittle uses unit diagonals for the lower triangle
template<typename Ty>
void Doolittle(std::vector<Ty>& S, std::vector<Ty>& D) {
	size_t d = size_t(std::sqrt(S.size()));
	assert(S.size() == d*d);
	assert(D.size() == d*d);
	for (size_t k = 0; k < d; ++k) {
		for (size_t j = k; j < d; ++j) {
			Ty sum = 0.;
			for (size_t p = 0; p < k; ++p) sum += D[k*d + p] * D[p*d + j];
			D[k*d + j] = (S[k*d + j] - sum); // not dividing by diagonals
		}
		for (size_t i = k + 1; i < d; ++i) {
			Ty sum = 0.;
			for (size_t p = 0; p < k; ++p) sum += D[i*d + p] * D[p*d + k];
			D[i*d + k] = (S[i*d + k] - sum) / D[k*d + k];
		}
	}
}
// SolveDoolittle takes an LU decomposition, LU, and a right hand side vector, b, and produces a result, x.
template<typename Ty>
void SolveDoolittle(const std::vector<Ty>& LU, const std::vector<Ty>& b, std::vector<Ty>& x) {
	int d = (int)b.size();
	std::vector<Ty> y(d);
	for (int i = 0; i < d; ++i) {
		Ty sum = 0.;
		for (int k = 0; k < i; ++k)sum += LU[i*d + k] * y[k];
		y[i] = (b[i] - sum); // not dividing by diagonals
	}
	for (int i = d - 1; i >= 0; --i) {
		Ty sum = 0.;
		for (int k = i + 1; k < d; ++k)sum += LU[i*d + k] * x[k];
		x[i] = (y[i] - sum) / LU[i*d + i];
	}
}

#endif

void GenerateTestCase(int a, int b) {
	std::cout << std::setw(3) << a << "/" << std::setw(3) << b << (sw::hprblas::isRepresentable(a, b) ? " is    " : " is not") << " representable " << (a / double(b)) << std::endl;
}

void EnumerateTestCases() {
	for (int i = 0; i < 30; i += 3) {
		for (int j = 0; j < 70; j += 7) {
			GenerateTestCase(i, j);
		}
	}
}

template<size_t nbits, size_t es, size_t capacity = 10>
void ComparePositDecompositions(mtl::mat::dense2D< sw::unum::posit<nbits, es> >& A, mtl::vec::dense_vector< sw::unum::posit<nbits, es> >& x, mtl::vec::dense_vector< sw::unum::posit<nbits, es> >& b) {
	using namespace sw::hprblas;
	size_t d = A.num_cols();
	assert(A.num_rows() == d);
	mtl::mat::dense2D< sw::unum::posit<nbits, es> > LU(d*d);

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
	size_t d = A.num_cols();
	assert(A.num_rows() == d);
	Matrix LU(d*d);

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
	int d = 5;

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
		mtl::mat::dense2D<IEEEType> Aieee(d*d);
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
		mtl::vec::dense_vector<IEEEType> bieee(d);
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
		mtl::mat::dense2D<PositType> Aposit(d*d);
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
		mtl::vec::dense_vector<PositType> bposit(d);
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
