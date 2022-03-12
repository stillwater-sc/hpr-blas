// performance.cpp: performance measurement of HPRBLAS algorithms
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

// Boost arbitrary precision ints
#include <boost/multiprecision/cpp_int.hpp>
#include <hprblas>

#include <universal/number/integer/integer.hpp>
#include <universal/number/cfloat/cfloat.hpp>
#include <universal/verification/performance_runner.hpp>

// measure performance of arithmetic operators
void TestUniversalIntegerOperatorPerformance() {
	using namespace sw::universal;
	std::cout << "\nUniversal Integer Arithmetic operator performance\n";

	uint64_t NR_OPS = 1000000;

	PerformanceRunner("integer<8>    add/subtract  ", AdditionSubtractionWorkload< sw::universal::integer<8> >, NR_OPS);
	PerformanceRunner("integer<16>   add/subtract  ", AdditionSubtractionWorkload< sw::universal::integer<16> >, NR_OPS);
	PerformanceRunner("integer<32>   add/subtract  ", AdditionSubtractionWorkload< sw::universal::integer<32> >, NR_OPS);
	PerformanceRunner("integer<64>   add/subtract  ", AdditionSubtractionWorkload< sw::universal::integer<64> >, NR_OPS);
	PerformanceRunner("integer<128>  add/subtract  ", AdditionSubtractionWorkload< sw::universal::integer<128> >, NR_OPS / 2);
	PerformanceRunner("integer<256>  add/subtract  ", AdditionSubtractionWorkload< sw::universal::integer<256> >, NR_OPS / 4);
	PerformanceRunner("integer<512>  add/subtract  ", AdditionSubtractionWorkload< sw::universal::integer<512> >, NR_OPS / 8);
	PerformanceRunner("integer<1024> add/subtract  ", AdditionSubtractionWorkload< sw::universal::integer<1024> >, NR_OPS / 16);

	NR_OPS = 1024 * 32;
	PerformanceRunner("integer<8>    division      ", DivisionWorkload< sw::universal::integer<8> >, NR_OPS);
	PerformanceRunner("integer<16>   division      ", DivisionWorkload< sw::universal::integer<16> >, NR_OPS);
	PerformanceRunner("integer<32>   division      ", DivisionWorkload< sw::universal::integer<32> >, NR_OPS);
	PerformanceRunner("integer<64>   division      ", DivisionWorkload< sw::universal::integer<64> >, NR_OPS / 2);
	PerformanceRunner("integer<128>  division      ", DivisionWorkload< sw::universal::integer<128> >, NR_OPS / 4);
	PerformanceRunner("integer<512>  division      ", DivisionWorkload< sw::universal::integer<512> >, NR_OPS / 8);
	PerformanceRunner("integer<1024> division      ", DivisionWorkload< sw::universal::integer<1024> >, NR_OPS / 16);

	NR_OPS = 1024 * 32;
	PerformanceRunner("integer<8>    multiplication", MultiplicationWorkload< sw::universal::integer<8> >, NR_OPS);
	PerformanceRunner("integer<16>   multiplication", MultiplicationWorkload< sw::universal::integer<16> >, NR_OPS);
	PerformanceRunner("integer<32>   multiplication", MultiplicationWorkload< sw::universal::integer<32> >, NR_OPS / 2);
	PerformanceRunner("integer<64>   multiplication", MultiplicationWorkload< sw::universal::integer<64> >, NR_OPS / 4);
	PerformanceRunner("integer<128>  multiplication", MultiplicationWorkload< sw::universal::integer<128> >, NR_OPS / 8);
	PerformanceRunner("integer<512>  multiplication", MultiplicationWorkload< sw::universal::integer<512> >, NR_OPS / 16);
	PerformanceRunner("integer<1024> multiplication", MultiplicationWorkload< sw::universal::integer<1024> >, NR_OPS / 32);
}

void TestBoostIntegerOperatorPerformance() {
	using namespace sw::universal;
	using namespace boost::multiprecision;
	std::cout << "\nUniversal Integer Arithmetic operator performance\n";


	uint64_t NR_OPS = 1000000;

//	PerformanceRunner("integer<8>    add/subtract  ", AdditionSubtractionWorkload< int8_t >, NR_OPS);
//	PerformanceRunner("integer<16>   add/subtract  ", AdditionSubtractionWorkload< int16_t>, NR_OPS);
//	PerformanceRunner("integer<32>   add/subtract  ", AdditionSubtractionWorkload< int32_t >, NR_OPS);
//	PerformanceRunner("integer<64>   add/subtract  ", AdditionSubtractionWorkload< int64_t >, NR_OPS);
	PerformanceRunner("integer<128>  add/subtract  ", AdditionSubtractionWorkload< int128_t >, NR_OPS / 2);
	PerformanceRunner("integer<256>  add/subtract  ", AdditionSubtractionWorkload< int256_t >, NR_OPS / 4);
	PerformanceRunner("integer<512>  add/subtract  ", AdditionSubtractionWorkload< int512_t >, NR_OPS / 8);
	PerformanceRunner("integer<1024> add/subtract  ", AdditionSubtractionWorkload< int1024_t >, NR_OPS / 16);

	NR_OPS = 1024 * 32;
//	PerformanceRunner("integer<8>    division      ", DivisionWorkload< int8_t>, NR_OPS);
//	PerformanceRunner("integer<16>   division      ", DivisionWorkload< int16_t >, NR_OPS);
//	PerformanceRunner("integer<32>   division      ", DivisionWorkload< int32_t >, NR_OPS);
//	PerformanceRunner("integer<64>   division      ", DivisionWorkload< int64_t >, NR_OPS / 2);
	PerformanceRunner("integer<128>  division      ", DivisionWorkload< int128_t >, NR_OPS / 4);
	PerformanceRunner("integer<512>  division      ", DivisionWorkload< int512_t >, NR_OPS / 8);
	PerformanceRunner("integer<1024> division      ", DivisionWorkload< int1024_t >, NR_OPS / 16);

	NR_OPS = 1024 * 32;
//	PerformanceRunner("integer<8>    multiplication", MultiplicationWorkload< int8_t >, NR_OPS);
//	PerformanceRunner("integer<16>   multiplication", MultiplicationWorkload< int16_t >, NR_OPS);
//	PerformanceRunner("integer<32>   multiplication", MultiplicationWorkload< int32_t >, NR_OPS / 2);
//	PerformanceRunner("integer<64>   multiplication", MultiplicationWorkload< int64_t >, NR_OPS / 4);
	PerformanceRunner("integer<128>  multiplication", MultiplicationWorkload< int128_t >, NR_OPS / 8);
	PerformanceRunner("integer<512>  multiplication", MultiplicationWorkload< int512_t >, NR_OPS / 16);
	PerformanceRunner("integer<1024> multiplication", MultiplicationWorkload< int1024_t >, NR_OPS / 32);
}

int main(int argc, char** argv)
try {
	using namespace std;
	using namespace mtl;
	using namespace sw::universal;
	using namespace sw::hprblas;
	using namespace boost::multiprecision;

	int nrOfFailedTestCases = 0;

#ifdef NOW
	cout << "Size of int128_t  : " << typeid(integer<128>).name() << " = " << sizeof(integer<128>) << endl;
	cout << "Size of int1024_t : " << typeid(integer<1024>).name() << " = " << sizeof(integer<1024>) << endl;

	// this type is not implemented
//	using spint32_t = number<cpp_int_backend<32, 32, signed_packed, unchecked, void> >;
//	cout << "Size of spint32_t : " << typeid(spint32_t).name() << " = " << sizeof(spint32_t) << endl;

	using mpint32_t = number<cpp_int_backend<32, 32, signed_magnitude, unchecked, void> >;
	cout << "Size of mpint32_t : " << typeid(mpint32_t).name() << " = " << sizeof(mpint32_t) << endl;
//	mpint32_t a; a = 10;

	using mpint64_t = number<cpp_int_backend<64, 64, signed_magnitude, unchecked, void> >;
	cout << "Size of mpint64_t : " << typeid(mpint64_t).name() << " = " << sizeof(mpint64_t) << endl;

	using checked_int64_t = number<cpp_int_backend<64, 64, signed_magnitude, checked, void> >;
	cout << "Size of chint64_t : " << typeid(checked_int64_t).name() << " = " << sizeof(checked_int64_t) << endl;

	cout << "Size of int128_t  : " << typeid(int128_t).name() << " = " << sizeof(int128_t) << endl;
	cout << "Size of int1024_t : " << typeid(int1024_t).name() << " = " << sizeof(int1024_t) << endl;

	{
		int1024_t a, b, c;
		a = int1024_t(1.0e20);
		b = int1024_t(2.0e20);
		c = a + b;
		cout << a << endl;
		cout << b << endl;
		cout << c << endl;
	}
#endif // NOW

	uint64_t NR_OPS = 1000000;
	PerformanceRunner("cfloat<128>   multiplication", MultiplicationWorkload< sw::universal::cfloat<128,15,uint32_t,true, false, false> >, NR_OPS / 8);
	PerformanceRunner("integer<128>  multiplication", MultiplicationWorkload< sw::universal::integer<128> >, NR_OPS / 8);
//	PerformanceRunner("integer<512>  multiplication", MultiplicationWorkload< sw::universal::integer<512> >, NR_OPS / 16);
//	PerformanceRunner("integer<1024> multiplication", MultiplicationWorkload< sw::universal::integer<1024> >, NR_OPS / 32);

//	TestUniversalIntegerOperatorPerformance();
//	TestBoostIntegerOperatorPerformance();

	return (nrOfFailedTestCases > 0 ? EXIT_FAILURE : EXIT_SUCCESS);
}
catch (char const* msg) {
	std::cerr << msg << std::endl;
	return EXIT_FAILURE;
}
catch (const sw::universal::universal_arithmetic_exception& err) {
	std::cerr << "Uncaught universal arithmetic exception: " << err.what() << std::endl;
	return EXIT_FAILURE;
}
catch (const sw::universal::universal_internal_exception& err) {
	std::cerr << "Uncaught universal internal exception: " << err.what() << std::endl;
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
