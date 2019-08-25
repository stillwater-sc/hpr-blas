// hprblas_c_api.cpp: C++ shim to deliver a C API for HPRBLAS
//
// Copyright (C) 2017-2019 Stillwater Supercomputing, Inc.
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.

#include <hprblas.h>

#ifdef __cplusplus
extern "C" {
#endif
	// fused matrix matrix multiply using posit8_t type
	void posit8_fmm(posit8_t* C, unsigned M, unsigned N, unsigned K, posit8_t* A, posit8_t* B) {
	}
	// fused matrix matrix multiply using posit16 type
	void posit16_fmm(posit16_t* C, unsigned M, unsigned N, unsigned K, posit16_t* A, posit16_t* B) {
	}
	// fused matrix matrix multiply using posit32 type
	void posit32_fmm(posit32_t* C, unsigned M, unsigned N, unsigned K, posit32_t* A, posit32_t* B) {
	}
	// fused matrix matrix multiply using posit64 type
	void posit64_fmm(posit64_t* C, unsigned M, unsigned N, unsigned K, posit64_t* A, posit64_t* B) {
	}
	// fused matrix matrix multiply using posit128 type
	void posit128_fmm(posit128_t* C, unsigned M, unsigned N, unsigned K, posit128_t* A, posit128_t* B) {
	}
	// fused matrix matrix multiply using posit256 type
	void posit256_fmm(posit256_t* C, unsigned M, unsigned N, unsigned K, posit256_t* A, posit256_t* B) {
	}
#ifdef __cplusplus
}
#endif