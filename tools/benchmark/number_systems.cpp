// number_systems.cpp: performance comparison between different number systems
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.
#include <boost/multiprecision/cpp_bin_float.hpp>
//#define POSIT_FAST_SPECIALIZATION 1
#define POSIT_FAST_POSIT_16_1 1
#define POSIT_FAST_POSIT_32_2 1
#include <universal/number/posit/posit>
#include <universal/performance/number_system.hpp>

// define a true 256-bit IEEE floating point type
constexpr size_t bits_in_octand = 113 + 128;
using cpp_bin_float_octand = boost::multiprecision::number<boost::multiprecision::backends::cpp_bin_float<bits_in_octand, boost::multiprecision::backends::digit_base_2, void, boost::int16_t, -16382, 16383>, boost::multiprecision::expression_template_option::et_off>;
// define the floating point types (single, double, quad, octand)
using sp = boost::multiprecision::cpp_bin_float_single;
using dp = boost::multiprecision::cpp_bin_float_double;
using qp = boost::multiprecision::cpp_bin_float_quad;
using op = cpp_bin_float_octand;

int main(int argc, char** argv)
try {
	using namespace std;
	using namespace sw::universal;

	cout << "Arithmetic performance comparison" << endl;

	posit<16,1> p16;
	float f;
	sp boostsp;
	posit<32, 2> p32;

	OperatorPerformance report;
	GeneratePerformanceReport(f, report);
	cout << ReportPerformance(f, report) << endl;
//	GeneratePerformanceReport(boostsp, report);
//	cout << ReportPerformance(boostsp, report) << endl;
	GeneratePerformanceReport(p16, report);
	cout << ReportPerformance(p16, report) << endl;
	GeneratePerformanceReport(p32, report);
	cout << ReportPerformance(p32, report) << endl;

	return EXIT_SUCCESS;
}
catch (char const* msg) {
	std::cerr << msg << std::endl;
	return EXIT_FAILURE;
}
catch (const sw::universal::posit_arithmetic_exception& err) {
	std::cerr << "Uncaught posit arithmetic exception: " << err.what() << std::endl;
	return EXIT_FAILURE;
}
catch (const sw::universal::quire_exception& err) {
	std::cerr << "Uncaught quire exception: " << err.what() << std::endl;
	return EXIT_FAILURE;
}
catch (const sw::universal::posit_internal_exception& err) {
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
