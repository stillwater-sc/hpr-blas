// two_sum.cpp: TwoSum evaluation of posit number systems
//
// Copyright (C) 2017-2019 Stillwater Supercomputing, Inc.
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.
#define ALIASING_ALLOWED
#include "common.hpp"
#include "../tests/test_helpers.hpp"
#include "../tests/posit_test_helpers.hpp"

/*
TwoSum denotes an algorithm introduced by Knuth in "The Art of Computer Programming", vol 2, Seminumerical Algorithms.

Given two floating point values a and b, generate a rounded sum s and a remainder r, such that
s = RoundToNearest(a + b), and
a + b = s + r

*/

// enumerate all addition cases for a posit configuration: is within 10sec till about nbits = 14
template<size_t nbits, size_t es>
int ValidateTwoSum(std::string tag, bool bReportIndividualTestCases) {
	using namespace std;
	const size_t NR_POSITS = (size_t(1) << nbits);
	int nrOfFailedTests = 0;
	using Posit = sw::unum::posit<nbits, es>;
	Posit pa, pb, psum, premainder, pref;
	pair<Posit, Posit> s_and_r;
	double da, db;
	for (size_t i = 0; i < NR_POSITS; i++) {
		pa.set_raw_bits(i);
		da = double(pa);
		for (size_t j = 0; j < NR_POSITS; j++) {
			pb.set_raw_bits(j);
			db = double(pb);

			s_and_r = sw::unum::twoSum(pa, pb);
			psum = s_and_r.first;
			premainder = s_and_r.second;
			pref = psum + premainder;
			psum = pa + pb;

			if (psum != pref) {
				nrOfFailedTests++;
				if (bReportIndividualTestCases)	ReportBinaryArithmeticError("FAIL", "+", pa, pb, pref, psum);
			}
			else {
				//if (bReportIndividualTestCases) ReportBinaryArithmeticSuccess("PASS", "+", pa, pb, pref, psum);
			}
		}
	}

	return nrOfFailedTests;
}

int main(int argc, char** argv)
try {
	using namespace std;
	using namespace sw::unum;

	constexpr size_t nbits = 32;
	constexpr size_t es = 2;
	using Posit = posit<nbits,es>;

	// print detailed bit-level computational intermediate results
	bool verbose = false;

	int nrOfFailedTestCases = 0;
	bool bReportIndividualTestCases = true;
	std::string tag = "TwoSum failed: ";

	// preserve the existing ostream precision
	auto precision = cout.precision();
	cout << setprecision(12);

	cout << "Posit TwoSum validation" << endl;

//	nrOfFailedTestCases += ReportTestResult(ValidateTwoSum<2, 0>(tag, bReportIndividualTestCases), "posit<2,0>", "twoSum");

//	nrOfFailedTestCases += ReportTestResult(ValidateTwoSum<3, 0>(tag, bReportIndividualTestCases), "posit<3,0>", "twoSum");
//	nrOfFailedTestCases += ReportTestResult(ValidateTwoSum<3, 1>(tag, bReportIndividualTestCases), "posit<3,1>", "twoSum");

	nrOfFailedTestCases += ReportTestResult(ValidateTwoSum<4, 0>(tag, bReportIndividualTestCases), "posit<4,0>", "twoSum");
	nrOfFailedTestCases += ReportTestResult(ValidateTwoSum<4, 1>(tag, bReportIndividualTestCases), "posit<4,1>", "twoSum");
	nrOfFailedTestCases += ReportTestResult(ValidateTwoSum<4, 2>(tag, bReportIndividualTestCases), "posit<4,2>", "twoSum");

	/*
	nrOfFailedTestCases += ReportTestResult(ValidateTwoSum<5, 0>(tag, bReportIndividualTestCases), "posit<5,0>", "twoSum");
	nrOfFailedTestCases += ReportTestResult(ValidateTwoSum<5, 1>(tag, bReportIndividualTestCases), "posit<5,1>", "twoSum");
	nrOfFailedTestCases += ReportTestResult(ValidateTwoSum<5, 2>(tag, bReportIndividualTestCases), "posit<5,2>", "twoSum");
	nrOfFailedTestCases += ReportTestResult(ValidateTwoSum<5, 3>(tag, bReportIndividualTestCases), "posit<5,3>", "twoSum");

	nrOfFailedTestCases += ReportTestResult(ValidateTwoSum<6, 0>(tag, bReportIndividualTestCases), "posit<6,0>", "twoSum");
	nrOfFailedTestCases += ReportTestResult(ValidateTwoSum<6, 1>(tag, bReportIndividualTestCases), "posit<6,1>", "twoSum");
	nrOfFailedTestCases += ReportTestResult(ValidateTwoSum<6, 2>(tag, bReportIndividualTestCases), "posit<6,2>", "twoSum");
	nrOfFailedTestCases += ReportTestResult(ValidateTwoSum<6, 3>(tag, bReportIndividualTestCases), "posit<6,3>", "twoSum");
	nrOfFailedTestCases += ReportTestResult(ValidateTwoSum<6, 4>(tag, bReportIndividualTestCases), "posit<6,4>", "twoSum");
	*/

	// restore the previous ostream precision
	cout << setprecision(precision);

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
