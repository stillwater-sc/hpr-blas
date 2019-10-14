// matmul.c: example C program using a fused matrix matrix multiply
//
// Copyright (C) 2017-2019 Stillwater Supercomputing, Inc.
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.

#include <hprblas.h>

typedef struct {
	unsigned rows;
	unsigned cols;
	unsigned sizeOfElement;
	void *data;
} Matrix;

Matrix* newMatrix(unsigned m, unsigned n, unsigned sizeOfElement) {
	Matrix *M = malloc(sizeof(Matrix));
	M->rows = m;
	M->cols = n;
	M->sizeOfElement = sizeOfElement;
	M->data = malloc(m * n * sizeOfElement);
	return M;
}

void clearMatrix(unsigned m, unsigned n, unsigned sizeOfElement, void* M) {
	for (unsigned i = 0; i < m; ++i) {

	}
}

int main(int argc, char* argv[]) {
	Matrix *A, *B, *C;

	int N = 5;
	A = newMatrix(5, 5, sizeof(double));
	B = newMatrix(5, 5, sizeof(double));
	C = newMatrix(5, 5, sizeof(double));

	// create input matrices


	// free matrices

}