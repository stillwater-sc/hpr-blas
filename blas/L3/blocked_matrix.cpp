// blocked_matrix.cpp: examples of using blocked matrix operations
//
// Copyright (C) 2017-2019 Stillwater Supercomputing, Inc.
//
// This file is part of the HPR-BLAS project, which is released under an MIT Open Source license.

// warning C4996: 'std::copy::_Unchecked_iterators::_Deprecate': Call to 'std::copy' with parameters that may be unsafe - this call relies on the caller to check that the passed values are correct.
// \mtl4\boost/numeric/mtl/operation/update.hpp(159): warning C4244: 'argument': conversion from 'const double' to 'float', possible loss of data
#pragma warning( disable : 4996 4244)
#include "common.hpp"
#include <hprblas>
// utilities to generate and print vectors and matrices
#include "utils/matvec.hpp"

template <typename Matrix>
void fill_and_print(Matrix& A, char name)
{
	// Set values in traditional way
	A = 1.2, 3.4,
		5.6, 7.8;

	// Just print them
	std::cout << name << " is \n" << A << "\n";
}

// ValidateBlocking iterates over block sizes to test the blocking algorithm
template<typename Scalar>
int ValidateBlocking(const mtl::mat::dense2D<Scalar> A, const mtl::mat::dense2D<Scalar>& B, const mtl::mat::dense2D<Scalar>& Reference) {
	unsigned nrOfTestFailures = 0;

	// validate preconditions for the function
	assert(mtl::mat::num_rows(A) == mtl::mat::num_cols(B));
	assert(mtl::mat::num_rows(A) == mtl::mat::num_cols(A));
	assert(mtl::mat::num_rows(B) == mtl::mat::num_cols(B));

	unsigned N = unsigned(mtl::mat::num_rows(A));

	// set up the validation iteration
	mtl::mat::dense2D<Scalar> C(N, N);
	for (unsigned blockSize = 2; blockSize < N; ++blockSize) {
		//C = sw::hprblas::bmm(A, B, blockSize);
		C = sw::hprblas::bfmm(A, B, blockSize);

		if (!sw::hprblas::isEqual(Reference, C)) {
			std::stringstream str;
			str << "Result matrix with block size " << blockSize;
			sw::hprblas::printMatrix(std::cout, str.str(), C);
			++nrOfTestFailures;
		}
		else {
			std::cout << "block size " << blockSize << " PASSED\n";
		}
	}
	return nrOfTestFailures;
}

#define MANUAL 0
#define STRESS 0

int main(int argc, char** argv)
try {
	using namespace std;
	using namespace sw::unum;
	using namespace mtl;
	using namespace sw::hprblas;
	
	constexpr size_t nbits = 8;
	constexpr size_t es = 0;
	using Scalar = posit<nbits, es>;
	using Quire = quire<nbits, es>;

	unsigned nrOfFailedTestCases = 0;
	unsigned blockSize = 3;
	unsigned N = 10;

#if MANUAL
	dense2D<Scalar> A(N,N), B(N,N);
	hessian_setup(A, 1);
//	A = 1;  // set to I
	B = 1;  // set to I

	dense2D<Scalar> subBlock(blockSize, blockSize);
	dense2D<Quire> partial(blockSize, blockSize);

	printSubMatrix(std::cout, "A(0,1)", A, 0, 1, blockSize, blockSize);
	printSubMatrix(std::cout, "B(1,1)", B, 1, 1, blockSize, blockSize);

	subBlock = Scalar(0);
	subBlockMM(subBlock, A, 0, 1, B, 1, 1);
	std::cout << subBlock << std::endl;

	return 0;

	subBlock = Scalar(0);
	for (unsigned bk = 0; bk < 2; ++bk) { // block iterator
		subBlockMM(subBlock, A, 0, bk, B, bk, 1);
		std::cout << subBlock << std::endl;
	}

#else
	dense2D<Scalar> A(N, N), B(N, N);
	hessian_setup(A, 1);
	B = 1;  // set to I

	mtl::mat::dense2D<Scalar> Reference(N, N);
	Reference = Scalar(1); // reference Identity matrix
	mtl::hessian_setup(Reference, 1);
	sw::hprblas::printMatrix(std::cout, "Reference Matrix", Reference);

	nrOfFailedTestCases = ValidateBlocking<Scalar>(A, B, Reference);
#endif // MANUAL

	if (nrOfFailedTestCases) cout << "FAIL" << endl; else cout << "PASS" << endl;

	return (nrOfFailedTestCases > 0 ? EXIT_FAILURE : EXIT_SUCCESS);
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
