// fused_dot.cpp example program showing a fused-dot product for error free linear algebra
//
// Copyright (C) 2017-2021 Stillwater Supercomputing, Inc.
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.
#include <chrono>
#include <hprblas>
// utilities to generate and print vectors and matrices
#include "utils/matvec.hpp"

constexpr double pi = 3.14159265358979323846;  // best practice for C++

/*
 We want to show the benefits of reproducibility provided by the FDP.

 When you have concurrent execution through blocked LA operators, you need control over the
 small numbers. There are three ways of doing it:
 1- sort the numbers: very expensive
 2- Kahn two_sum accumulation
 3- fused-dot product using a super-accumulator (Ulrich Kulisch)

 We can construct a tough test case, where we have segments of decreasing numbers along the accumulation path.
 That will allow us to construct specific cancellation cases that yield nice predictable answers.

 We would also like to collect real-world cases of dot products that go wrong when not accounting for
 catastrophic cancellation possibilities.

 */

template<typename Vector>
void FillDescending(Vector& vec, typename Vector::value_type start) {
	using namespace mtl;
	size_t n = mtl::size(vec);
	//std::cout << "Nr of samples to generate is " << n << std::endl;
	for (size_t i = 0; i < n; ++i) {
		vec[i] = start;
		start *= 0.5;
	}
}

template<typename Vector>
void FillAscending(Vector& vec, typename Vector::value_type start) {
	using namespace mtl;
	size_t n = mtl::size(vec);
	//std::cout << "Nr of samples to generate is " << n << std::endl;
	size_t r = n;
	for (size_t i = 0; i < n; ++i) {
		vec[--r] = start;
		start *= 0.5;
	}
}

int main(int argc, char** argv)
try {
	using namespace std::chrono;
	using namespace mtl;
	using namespace sw::universal;
	using namespace sw::hprblas;

	constexpr size_t nbits = 16;
	constexpr size_t es = 1;
	constexpr size_t vecSize = 32;

	int nrOfFailedTestCases = 0;

	{
		using Scalar = posit<nbits, es>;
		Scalar pmaxpos(SpecificValue::maxpos);
		std::vector<Scalar> px(vecSize), py(vecSize);

		FillDescending(px, 0.25*pmaxpos);
		FillAscending(py, 0.25*pmaxpos);

		printVector(std::cout, "px", px);
		printVector(std::cout, "px", py);

		steady_clock::time_point t1 = steady_clock::now();
		Scalar presult = sw::hprblas::fdp(px, py);
		steady_clock::time_point t2 = steady_clock::now();
		double ops = vecSize * 2.0; // dot product is vecSize products and vecSize adds
		std::cout << "FDP product     is " << presult << '\n';
		duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
		double elapsed = time_span.count();
		std::cout << "It took " << elapsed << " seconds." << std::endl;
		std::cout << "Performance " << (uint32_t)(ops / (1000 * elapsed)) << " KOPS" << '\n';
	}

	{
		using Scalar = posit<2*nbits, es+1>;
		Scalar pmaxpos(SpecificValue::maxpos);
		std::vector<Scalar> px(vecSize), py(vecSize);

		FillDescending(px, 0.25*pmaxpos);
		FillAscending(py, 0.25*pmaxpos);

//		printVector(cout, "px", px);
//		printVector(cout, "px", py);

		steady_clock::time_point t1 = steady_clock::now();
		Scalar presult = sw::hprblas::dot(px, py);
		steady_clock::time_point t2 = steady_clock::now();
		double ops = vecSize * 2.0; // dot product is vecSize products and vecSize adds
		std::cout << "DOT product     is " << presult << '\n';
		duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
		double elapsed = time_span.count();
		std::cout << "It took " << elapsed << " seconds." << std::endl;
		std::cout << "Performance " << (uint32_t)(ops / (1000 * elapsed)) << " KOPS" << '\n';
	}

	{
		using Scalar = float;
		std::vector<Scalar> px(vecSize), py(vecSize);

		posit<32,2> pmaxpos(SpecificValue::maxpos);
		FillDescending(px, float(0.25*pmaxpos));
		FillAscending(py, float(0.25*pmaxpos));

		//		printVector(cout, "px", px);
		//		printVector(cout, "px", py);

		steady_clock::time_point t1 = steady_clock::now();
		Scalar presult = sw::hprblas::dot(px, py);
		steady_clock::time_point t2 = steady_clock::now();
		double ops = vecSize * 2.0; // dot product is vecSize products and vecSize adds
		std::cout << "DOT product     is " << presult << '\n';
		duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
		double elapsed = time_span.count();
		std::cout << "It took " << elapsed << " seconds." << std::endl;
		std::cout << "Performance " << (uint32_t)(ops / (1000 * elapsed)) << " KOPS" << '\n';
	}

	return (nrOfFailedTestCases > 0 ? EXIT_FAILURE : EXIT_SUCCESS);
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
