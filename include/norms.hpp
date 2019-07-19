#pragma once
// norms.hpp : vector and matrix norms
//
// Copyright (C) 2017-2019 Stillwater Supercomputing, Inc.
//
// This file is part of the HPR-BLAS project, which is released under an MIT Open Source license.

/*
a general vector norm | x | , sometimes written with a double bar as || x ||, 
is a nonnegative norm defined such that

1. | x | > 0 when x != 0 and | x |= 0 iff x = 0.

2. | kx |= | k || x | for any scalar k.

3. | x + y |<= |x| + |y|
*/

namespace sw {
namespace hprblas {

// L1-norm = sum of absolute value of each vector element
template<typename Vector>
typename Vector::value_type l1_norm(const Vector& v) {
	using Scalar = typename Vector::value_type;

	Scalar l1 = Scalar(0);
	for (unsigned i = 0; i < size(v); ++i) {
		l1 += abs(v[i]);
	}
	return l1;
}

// L2-norm = Manhattan distance
template<typename Vector>
typename Vector::value_type l2_norm(const Vector& v) {
	using Scalar = typename Vector::value_type;

	Scalar l2 = Scalar(0);
	for (unsigned i = 0; i < size(v); ++i) {
		l2 += v[i] * v[i];
	}
	return sqrt(l2);
}

// Linfinity-norm = max of the absolute value of each vector element
template<typename Vector>
typename Vector::value_type linf_norm(const Vector& v) {
	using Scalar = typename Vector::value_type;

	Scalar linf = Scalar(0);
	for (unsigned i = 0; i < size(v); ++i) {
		linf = abs(v[i]) > linf ? abs(v[i]) : linf;
	}
	return linf;
}

// calculate an integer power function base^positive_int
template<typename Scalar>
Scalar integer_power(Scalar base, int exponent) {
	if (exponent < 0) {
		base = Scalar(1) / base;
		exponent = -exponent;
	}
	if (exponent == 0) return Scalar(1);
	Scalar power = Scalar(1);
	while (exponent > 1) {
		if (exponent & 0x1) {
			power = base * power;
			base *= base;
			exponent = (exponent - 1) / 2;
		}
		else {
			base *= base;
			exponent /= 2;
		}
	}
	return base * power;
}

template<typename Scalar>
Scalar error_volume(Scalar Linf, unsigned dimensionality, bool measuredInULP = false) {
	Scalar ulp = std::numeric_limits<Scalar>::epsilon();
	if (measuredInULP) {
		return integer_power(Linf / ulp, dimensionality);
	}
	return integer_power(Linf, dimensionality);
}

} // namespace hprblas
} // namespace sw
