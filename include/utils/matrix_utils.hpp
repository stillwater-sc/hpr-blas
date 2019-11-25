#pragma once
// matrix_utils.hpp :  include file containing templated utilities to work with matrices
//
// Copyright (C) 2017-2019 Stillwater Supercomputing, Inc.
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.
#include <cstdint>
#include <random>
#include <algorithm>

namespace sw {
namespace hprblas {

// fill a dense matrix with random values between [lowerbound, upperbound]
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

	typedef typename mtl::Collection<Matrix>::value_type    value_type;
	typedef typename mtl::Collection<Matrix>::size_type     size_type;

	// inserters add to the elements, so we need to set the value to 0 before we begin
	A = 0.0;
	{ // extra block unfortunately needed for VS2013
		// Create inserter for matrix m
		mtl::mat::inserter<Matrix> ins(A, num_cols(A));

		// generate and insert random values in A
		for (size_type r = 0; r < num_rows(A); r++) {
			for (size_type c = 0; c < num_cols(A); c++) {
				ins[r][c] << value_type(dist(engine));
			}
		}
		// Destructor of ins sets final state of m
	}
}

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

// fill a dense MTL matrix with diagonally dominant random values between [lowerbound, upperbound]
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

// Random Orthogonal Matrices

/*
Standard methods for generating random orthogonal matrices with Haar distribution 
are based on the method of Heiberger (1978). With this method, an (n x n) matrix A is
first generated with entries xij ~ Normal(0; 1). Then a QR factorization (X = QR) is computed. 
This method provides a random Q with correct distribution.
*/
template<typename Matrix>
void uniform_random_orthogonal_Heiberger(Matrix& Q) {
	// Use random_device to generate a seed for Mersenne twister engine.
	std::random_device rd{};
	// Use Mersenne twister engine to generate pseudo-random numbers.
	std::mt19937 engine{ rd() };
	// "Filter" MT engine's output to generate pseudo-random double values,
	// **uniformly distributed** on the closed interval [lowerbound, upperbound].
	// (Note that the range is [inclusive, inclusive].)
	std::uniform_real_distribution<double> dist{ 0.0, 1.0 };
	// Pattern to generate pseudo-random number.
	// double rnd_value = dist(engine);

	typedef typename mtl::Collection<Matrix>::value_type    value_type;
	typedef typename mtl::Collection<Matrix>::size_type     size_type;
	Matrix A(mtl::mat::num_rows(Q), mtl::mat::num_cols(Q));
	Matrix R(mtl::mat::num_rows(Q), mtl::mat::num_cols(Q));

	// fill X with elements from Normal(0,1)
	// inserters add to the elements, so we need to set the value to 0 before we begin
	A = 0.0;
	{ // extra block unfortunately needed for VS2013
		// Create inserter for matrix A
		mtl::mat::inserter<Matrix> ins(A, num_cols(A));

		// insert random values in A
		for (size_type r = 0; r < num_rows(A); r++) {
			for (size_type c = 0; c < num_cols(A); c++) {
				ins[r][c] << value_type(dist(engine));
			}
		}
		// Destructor of ins sets final state of A
	}
	// QR method to generate orthonormal matrix Q
	mtl::mat::qr(A, Q, R);
	std::cout << A << std::endl;
	std::cout << Q << std::endl;
	std::cout << R << std::endl;
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
// m[i,j]= {  1 if A[i,j]>=0
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
void uniform_rand_orthogonal(Matrix& A) {

}

//
// fill a dense upper-triangular matrix with elements: [i][j] = 0.5 + (1 + j - i)*0.5
template <typename Matrix>
void fill_U(Matrix& A, double lowerbound = 0.0, double upperbound = 1.0)
{
	using namespace mtl;
	typedef typename Collection<Matrix>::value_type    value_type;
	typedef typename Collection<Matrix>::size_type     size_type;

	// inserters add to the elements, so we need to set the value to 0 before we begin
	A = 0.0;
	{ // extra block unfortunately needed for VS2013
		// Create inserter for matrix m
		mat::inserter<Matrix> ins(A, num_cols(A));

		// generate and insert random values in A
		for (size_type r = 0; r < num_rows(A); r++) {
			for (size_type c = r; c < num_cols(A); c++) {
				ins[r][c] << value_type(0.5) + value_type(1 + c - r)*value_type(0.5);
				ins[r][c] << value_type(1.0);
			}
		}
		// Destructor of ins sets final state of m
	}
}

// fill a dense lower-triangular matrix with elements: [i][j] = 0.5 + (1 + i - j)*0.5
template <typename Matrix>
void fill_L(Matrix& A, double lowerbound = 0.0, double upperbound = 1.0)
{
	using namespace mtl;
	typedef typename Collection<Matrix>::value_type    value_type;
	typedef typename Collection<Matrix>::size_type     size_type;

	// inserters add to the elements, so we need to set the value to 0 before we begin
	A = 0.0;
	{ // extra block unfortunately needed for VS2013
		// Create inserter for matrix m
		mat::inserter<Matrix> ins(A, num_cols(A));

		// generate and insert random values in A
		for (size_type r = 0; r < num_rows(A); r++) {
			for (size_type c = 0; c <= r; c++) {
				ins[r][c] << value_type(0.5) + value_type(1 + r - c)*value_type(0.5);
				ins[r][c] << value_type(1.0);

			}
		}
		// Destructor of ins sets final state of m
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
/// Hilbert matrices

// Greatest Common Divisor of two numbers, a and b
template<typename IntegerType>
IntegerType gcd(IntegerType a, IntegerType b)
{
	if (b == IntegerType(0))	return a;
	return gcd(b, a % b);
}

// Least Common Multiple of n numbers
template<typename IntegerType>
IntegerType findlcm(const std::vector<IntegerType>& v)
{
	IntegerType lcm = v[0];
	for (size_t i = 1; i < v.size(); i++) {
		lcm = (v[i] * lcm) / gcd(v[i], lcm);
	}
	return lcm;
}

// Generate the scaling factor of a Hilbert matrix so that its elements are representable
// that is, no infinite expensions of rationals, such as 1/3, 1/10, etc.
template<typename IntegerType>
IntegerType HilbertScalingFactor(IntegerType N) {
	std::vector<IntegerType> coef;
	for (IntegerType i = 2; i <= N; ++i) coef.push_back(i);
	for (IntegerType j = 2; j <= N; ++j) coef.push_back(N + j - IntegerType(1));
	return findlcm(coef);
}

template<typename Scalar>
size_t GenerateHilbertMatrix(mtl::mat::dense2D<Scalar>& m, bool bScale = true) {
	assert(m.num_rows() == m.num_cols());
	size_t N = m.num_rows();
	size_t lcm = HilbertScalingFactor(N); // always calculate the Least Common Multiplier
	Scalar scale = bScale ? lcm : Scalar(1);
	for (int i = 1; i <= N; ++i) {
		for (int j = 1; j <= N; ++j) {
			m[i - 1][j - 1] = scale / Scalar(i + j - 1);
		}
	}
	return lcm;
}

uint64_t factorial(uint64_t n) {
	return (n == 0 || n == 1) ? 1 : factorial(n - 1) * n;
}

// (n over k) = n! / k!(n-k)!
// this is computationally difficult due to the large numbers generated by the factorials
template<typename Scalar>
Scalar BinomialCoefficient_factorial(uint64_t n, uint64_t k) {
	Scalar numerator = Scalar(factorial(n));
	Scalar denominator = Scalar(factorial(k)*factorial(n - k));
	Scalar coef = numerator / denominator;
	//	std::cout << numerator << " / " << denominator << " = " << coef << std::endl;
	if (coef * denominator != numerator) std::cout << "FAIL: (" << n << " over " << k << ")  coef = " << coef << " ( " << numerator << "/" << denominator << ") error = " << numerator - (coef * denominator) << std::endl;
	return coef;
}

// (n over k) = (n-1 over k-1) + (n-1 over k)
template<typename Scalar>
Scalar BinomialCoefficient(Scalar n, Scalar k) {
	Scalar coef;
	if (k == 0 || k == n || n == 0) {
		coef = Scalar(1);
	}
	else {
		coef = BinomialCoefficient(n - 1, k - 1) + BinomialCoefficient(n - 1, k);
	}
	return coef;
}

template<typename Scalar>
void GenerateHilbertMatrixInverse(mtl::mat::dense2D<Scalar>& m, Scalar scale = Scalar(1.0)) {
	assert(m.num_rows() == m.num_cols());
	size_t N = m.num_rows();
	for (int i = 1; i <= N; ++i) {
		for (int j = 1; j <= N; ++j) {
			Scalar sign = ((i + j) % 2) ? Scalar(-1) : Scalar(1);
			Scalar factor1 = Scalar(i + j - 1);
			Scalar factor2 = Scalar(BinomialCoefficient<uint64_t>(N + i - 1, N - j));
			Scalar factor3 = Scalar(BinomialCoefficient<uint64_t>(N + j - 1, N - i));
			Scalar factor4 = Scalar(BinomialCoefficient<uint64_t>(i + j - 2, i - 1));
			m[i - 1][j - 1] = Scalar(sign * factor1 * factor2 * factor3 * factor4 * factor4);
			/* for tracing dynamic range failures
			std::cout << "element " << i << "," << j << std::endl;
			std::cout << "sign    " << sign << std::endl;
			std::cout << "factor1 " << factor1 << std::endl;
			std::cout << "factor2 " << factor2 << std::endl;
			std::cout << "factor3 " << factor3 << std::endl;
			std::cout << "factor4 " << factor4 << std::endl;
			*/
		}
	}
}

// isEqual compares to matrices
template<typename Matrix>
bool isEqual(const Matrix& lhs, const Matrix& rhs) {
	size_t r = lhs.num_rows();
	size_t c = lhs.num_cols();
	for (size_t i = 0; i < r; ++i) {
		for (size_t j = 0; j < c; ++j) {
			if (lhs[i][j] != rhs[i][j]) return false;
		}
	}
	return true;
}

} // namespace hprblas
} // namespace sw