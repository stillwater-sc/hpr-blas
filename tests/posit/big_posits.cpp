// big_posit.cpp: Functionality tests for big posits
//
// Copyright (C) 2017-2019 Stillwater Supercomputing, Inc.
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.

#include "common.hpp"

#include <boost/multiprecision/cpp_bin_float.hpp>
// enable posit arithmetic exceptions
#define POSIT_THROW_ARITHMETIC_EXCEPTION 1
// to capture all the possible bits, use 
#define POSIT_ROUNDING_ERROR_FREE_IO_FORMAT 1
#include <posit>
#include "../test_helpers.hpp"
#include "../posit_test_helpers.hpp"

/*
experiments with big posits
*/

// Sample the conversion space around 1 when number of fraction bits
// of a big posit is bigger than the number of bits in the input representation
template<typename Ty>
void Sample(Ty value) {
	using namespace std;
	using namespace sw::unum;

	cout << typeid(value).name() << endl;
	cout << posit<54, 3>(value) << " " << double(posit<54, 3>(value)) << endl;
	cout << posit<56, 3>(value) << " " << double(posit<56, 3>(value)) << endl;
	cout << posit<58, 3>(value) << " " << double(posit<58, 3>(value)) << endl;
	cout << posit<60, 3>(value) << " " << double(posit<60, 3>(value)) << endl;
	cout << posit<62, 3>(value) << " " << double(posit<62, 3>(value)) << endl;
	cout << posit<64, 3>(value) << " " << double(posit<64, 3>(value)) << endl;
	cout << posit<66, 3>(value) << " " << double(posit<66, 3>(value)) << endl;
	cout << posit<67, 3>(value) << " " << double(posit<67, 3>(value)) << endl;
	cout << posit<68, 3>(value) << " " << double(posit<68, 3>(value)) << endl;
	cout << posit<69, 3>(value) << " " << double(posit<69, 3>(value)) << endl;
	cout << posit<70, 3>(value) << " " << double(posit<70, 3>(value)) << endl;
	cout << posit<71, 3>(value) << " " << double(posit<71, 3>(value)) << endl;
	cout << posit<72, 3>(value) << " " << double(posit<72, 3>(value)) << endl;
	cout << posit<80, 3>(value) << " " << double(posit<80, 3>(value)) << endl;
}

int main(int argc, char** argv)
try {
	using namespace std;
	using namespace sw::unum;

	const size_t RND_TEST_CASES = 1000;

	const size_t nbits = 128;
	const size_t es = 4;

	int nrOfFailedTestCases = 0;
	bool bReportIndividualTestCases = false;
	std::string tag = " big posit conversion experiments";

	Sample<unsigned long>((unsigned long)1);
	Sample<unsigned long>((unsigned long)2);
	Sample<unsigned long long>((unsigned long long)1);
	Sample<unsigned long long>((unsigned long long)2);
	Sample<unsigned long long>((unsigned long long)4);

	return (nrOfFailedTestCases > 0 ? EXIT_FAILURE : EXIT_SUCCESS);
}
catch (char const* msg) {
	std::cerr << msg << std::endl;
	return EXIT_SUCCESS; //as we manually throwing the not supported yet it should not fall through the cracks     EXIT_FAILURE;
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
catch (const std::runtime_error& err) {
	std::cerr << "Uncaught runtime exception: " << err.what() << std::endl;
	return EXIT_FAILURE;
}
catch (...) {
	std::cerr << "Caught unknown exception" << std::endl;
	return EXIT_FAILURE;
}