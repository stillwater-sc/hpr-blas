// performance.cpp : performance benchmarking for abitrary fixed-precision cfloats
//
// Copyright (C) 2017-2021 Stillwater Supercomputing, Inc.
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.
#include <universal/utility/directives.hpp>
#include <chrono>
#include <boost/multiprecision/cpp_bin_float.hpp>
// configure the cfloat arithmetic class
#define cfloat_THROW_ARITHMETIC_EXCEPTION 0
#include <universal/number/cfloat/cfloat.hpp>
// is representable
#include <universal/functions/isrepresentable.hpp>
#include <universal/verification/test_status.hpp> // ReportTestResult
#include <universal/verification/performance_runner.hpp>

/*
   The goal of the arbitrary fixed-precision cfloats is to provide a constrained 
   linear floating-point type to explore the benefits of mixed-precision algorithms.
*/

template<typename cfloatConfiguration>
void CopyWorkload(size_t NR_OPS) {
	using namespace sw::universal;
	cfloatConfiguration a,b,c;

	bool bFail = false;
	size_t j = 0;
	for (size_t i = 0; i < NR_OPS; ++i,++j) {
		a.setbits(i);
		b = a;
		c.setbits(j);
		if (b.sign() != c.sign()) {
			bFail = true;
		}
	}
	if (bFail) std::cout << "COPY FAIL\n"; // just a quick double check that all went well
}


/// <summary>
/// measure performance of copying numbers around
/// </summary>
void TestCopyPerformance() {
	using namespace boost::multiprecision;
	using namespace sw::universal;
	std::cout << "copy performance comparisoni\n";

	uint64_t NR_OPS = 10000000;

	std::cout << "very large representations\n";
	PerformanceRunner("cfloat<80,11,uint64_t>   copy           ", CopyWorkload< sw::universal::cfloat<80, 11, uint64_t> >, NR_OPS);
	PerformanceRunner("cfloat<96,11,uint64_t>   copy           ", CopyWorkload< sw::universal::cfloat<96, 11, uint64_t> >, NR_OPS);
	PerformanceRunner("cfloat<128,11,uint64_t>  copy           ", CopyWorkload< sw::universal::cfloat<128, 11, uint64_t> >, NR_OPS);
	PerformanceRunner("cfloat<256,11,uint64_t>  copy           ", CopyWorkload< sw::universal::cfloat<256, 11, uint64_t> >, NR_OPS);
	PerformanceRunner("cfloat<512,11,uint64_t>  copy           ", CopyWorkload< sw::universal::cfloat<512, 11, uint64_t> >, NR_OPS);
}

template<typename Scalar>
void DecodeWorkload(size_t NR_OPS) {
	using namespace sw::universal;

	Scalar a{ 0 };
	size_t success{ 0 };
	bool first{ true };
	for (size_t i = 0; i < NR_OPS; ++i) {
		a.setbits(i);
		bool s{ false };
		blockbinary<a.es, typename Scalar::BlockType> e;
		blockbinary<a.fbits, typename Scalar::BlockType> f;
		sw::universal::decode(a, s, e, f);
		if ((i%2) == f.at(0)) {
			++success;
		}
		else {
			// this shouldn't happen, but found a bug this way with cfloat<64,11,uint64_t> as type
			if (first) {
				first = false;
				std::cout << typeid(a).name() << " :\n"
					<< to_binary(a,true) << "\n" 
					<< "sign    : " << (s ? "-1\n" : "+1\n") 
					<< "exponent: " << to_binary(e,true) << "\n" 
					<< "fraction: " << to_binary(f,true) << "\n";
			}
		}
	}
	if (success == 0) std::cout << "DECODE FAIL\n"; // just a quick double check that all went well
}


/// <summary>
/// measure performance of decode operator
/// NOTE: es is <= 11 due to limits of dynamic range of a 64-bit double
/// </summary>
void TestDecodePerformance() {
	using namespace sw::universal;
	std::cout << "classic floating-point cfloat decode operator performance\n";

	uint64_t NR_OPS = 100000;


	std::cout << "very large representations\n";
	PerformanceRunner("cfloat<64,11,uint64_t>   decode         ", DecodeWorkload< sw::universal::cfloat<64, 11, uint64_t> >, NR_OPS);
	PerformanceRunner("cfloat<80,11,uint64_t>   decode         ", DecodeWorkload< sw::universal::cfloat<80, 11, uint32_t> >, NR_OPS);
	PerformanceRunner("cfloat<96,11,uint64_t>   decode         ", DecodeWorkload< sw::universal::cfloat<96, 11, uint32_t> >, NR_OPS);
	PerformanceRunner("cfloat<128,15,uint64_t>  decode         ", DecodeWorkload< sw::universal::cfloat<128, 15, uint32_t> >, NR_OPS);
	PerformanceRunner("cfloat<256,15,uint64_t>  decode         ", DecodeWorkload< sw::universal::cfloat<256, 15, uint32_t> >, NR_OPS);
	PerformanceRunner("cfloat<512,15,uint64_t>  decode         ", DecodeWorkload< sw::universal::cfloat<512, 15, uint32_t> >, NR_OPS);
}

// measure performance of conversion operators
void TestConversionPerformance() {
	using namespace sw::universal;
	std::cout << "classic floating-point cfloat conversion performance\n";

//	uint64_t NR_OPS = 1000000;
}

// measure performance of arithmetic operators
void TestArithmeticOperatorPerformance() {
	using namespace sw::universal;
	std::cout << "classic floating-point cfloat arithmetic operator performance\n";

	uint64_t NR_OPS = 1000000;

	PerformanceRunner("cfloat<8,2,uint8_t>      add/subtract   ", AdditionSubtractionWorkload< sw::universal::cfloat<8,2,uint8_t> >, NR_OPS);
	PerformanceRunner("cfloat<16,5,uint16_t>    add/subtract   ", AdditionSubtractionWorkload< sw::universal::cfloat<16,5,uint16_t> >, NR_OPS);
	PerformanceRunner("cfloat<32,8,uint32_t>    add/subtract   ", AdditionSubtractionWorkload< sw::universal::cfloat<32,8,uint32_t> >, NR_OPS);
	PerformanceRunner("cfloat<64,11,uint64_t>   add/subtract   ", AdditionSubtractionWorkload< sw::universal::cfloat<64,11,uint64_t> >, NR_OPS);
	PerformanceRunner("cfloat<128,11,uint32_t>  add/subtract   ", AdditionSubtractionWorkload< sw::universal::cfloat<128,11, uint32_t> >, NR_OPS / 2);
	PerformanceRunner("cfloat<128,15,uint32_t>  add/subtract   ", AdditionSubtractionWorkload< sw::universal::cfloat<128,15, uint32_t> >, NR_OPS / 2);
	PerformanceRunner("cfloat<256,15,uint32_t   add/subtract   ", AdditionSubtractionWorkload< sw::universal::cfloat<256,15, uint32_t> >, NR_OPS / 4);
	PerformanceRunner("cfloat<512,15,uint32_t>  add/subtract   ", AdditionSubtractionWorkload< sw::universal::cfloat<512,15, uint32_t> >, NR_OPS / 8);
	PerformanceRunner("cfloat<1024,15,uint32_t> add/subtract   ", AdditionSubtractionWorkload< sw::universal::cfloat<1024,15, uint32_t> >, NR_OPS / 16);

	NR_OPS = 1024 * 32;
	PerformanceRunner("cfloat<8,2,uint16_t>     division       ", DivisionWorkload< sw::universal::cfloat<8,2,uint16_t> >, NR_OPS);
	PerformanceRunner("cfloat<16,5,uint16_t>    division       ", DivisionWorkload< sw::universal::cfloat<16,5,uint16_t> >, NR_OPS);
	PerformanceRunner("cfloat<32,8,uint32_t>    division       ", DivisionWorkload< sw::universal::cfloat<32,8,uint32_t> >, NR_OPS);
	PerformanceRunner("cfloat<64,11,uint64_t>   division       ", DivisionWorkload< sw::universal::cfloat<64,11,uint64_t> >, NR_OPS);
	PerformanceRunner("cfloat<128,15,uint32_t>  division       ", DivisionWorkload< sw::universal::cfloat<128,15, uint32_t> >, NR_OPS / 2);
	PerformanceRunner("cfloat<256,15,uint32_t   division       ", DivisionWorkload< sw::universal::cfloat<256,15, uint32_t> >, NR_OPS / 4);
	PerformanceRunner("cfloat<512,15,uint32_t>  division       ", DivisionWorkload< sw::universal::cfloat<512,15, uint32_t> >, NR_OPS / 8);
	PerformanceRunner("cfloat<1024,15,uint32_t> division       ", DivisionWorkload< sw::universal::cfloat<1024,15, uint32_t> >, NR_OPS / 16);

	// multiplication is the slowest operator

	NR_OPS = 1024 * 32;
	PerformanceRunner("cfloat<8,2,uint16_t>     multiplication ", MultiplicationWorkload< sw::universal::cfloat<8,2,uint16_t> >, NR_OPS);
	PerformanceRunner("cfloat<16,5,uint16_t>    multiplication ", MultiplicationWorkload< sw::universal::cfloat<16,5,uint16_t> >, NR_OPS);
	PerformanceRunner("cfloat<32,8,uint32_t>    multiplication ", MultiplicationWorkload< sw::universal::cfloat<32,8,uint32_t> >, NR_OPS);
	PerformanceRunner("cfloat<64,11,uint64_t>   multiplication ", MultiplicationWorkload< sw::universal::cfloat<64,11,uint64_t> >, NR_OPS);
	PerformanceRunner("cfloat<128,15,uint32_t>  multiplication ", MultiplicationWorkload< sw::universal::cfloat<128,15, uint32_t> >, NR_OPS / 2);
	PerformanceRunner("cfloat<256,15,uint32_t   multiplication ", MultiplicationWorkload< sw::universal::cfloat<256,15, uint32_t> >, NR_OPS / 4);
	PerformanceRunner("cfloat<512,15,uint32_t>  multiplication ", MultiplicationWorkload< sw::universal::cfloat<512,15, uint32_t> >, NR_OPS / 8);
	PerformanceRunner("cfloat<1024,15,uint32_t> multiplication ", MultiplicationWorkload< sw::universal::cfloat<1024,15, uint32_t> >, NR_OPS / 16);

}

// measure performance of arithmetic operators
void TestBoostMPArithmeticOperatorPerformance() {
	using namespace sw::universal;
	using namespace boost::multiprecision;
	std::cout << "boost::multiprecision arithmetic operator performance\n";

	using boostmp_8_2_uint8_t = number<backends::cpp_bin_float<5, backends::digit_base_2, void, std::int8_t, -0, 1>, et_off>;
	using boostmp_16_5_uint16_t = number<backends::cpp_bin_float<5, backends::digit_base_2, void, std::int16_t, -14, 15>, et_off>;
	using boostmp_32_8_uint32_t = number<backends::cpp_bin_float<24, backends::digit_base_2, void, std::int16_t, -126, 127>, et_off>;
	using boostmp_64_11_uint64_t = number<backends::cpp_bin_float<52, backends::digit_base_2, void, std::int16_t, -1022, 1023>, et_off>;
	using boostmp_128_11_uint32_t = number<backends::cpp_bin_float<116, backends::digit_base_2, void, std::int16_t, -1022, 1023>, et_off>;
	using boostmp_128_15_uint32_t = number<backends::cpp_bin_float<112, backends::digit_base_2, void, std::int16_t, -16382, 16385>, et_off>;
	using boostmp_256_15_uint32_t = number<backends::cpp_bin_float<240, backends::digit_base_2, void, std::int16_t, -16382, 16385>, et_off>;
	using boostmp_512_15_uint32_t = number<backends::cpp_bin_float<496, backends::digit_base_2, void, std::int16_t, -16382, 16385>, et_off>;
	using boostmp_1024_15_uint32_t = number<backends::cpp_bin_float<1008, backends::digit_base_2, void, std::int16_t, -16382, 16385>, et_off>;

	uint64_t NR_OPS = 1000000;

	PerformanceRunner("boostmp<8,2,uint8_t>      add/subtract   ", AdditionSubtractionWorkload< boostmp_8_2_uint8_t >, NR_OPS);
	PerformanceRunner("boostmp<16,5,uint16_t>    add/subtract   ", AdditionSubtractionWorkload< boostmp_16_5_uint16_t >, NR_OPS);
	PerformanceRunner("boostmp<32,8,uint32_t>    add/subtract   ", AdditionSubtractionWorkload< boostmp_32_8_uint32_t >, NR_OPS);
	PerformanceRunner("boostmp<64,11,uint64_t>   add/subtract   ", AdditionSubtractionWorkload< boostmp_64_11_uint64_t >, NR_OPS);
	PerformanceRunner("boostmp<128,11,uint32_t>  add/subtract   ", AdditionSubtractionWorkload< boostmp_128_11_uint32_t >, NR_OPS / 2);
	PerformanceRunner("boostmp<128,15,uint32_t>  add/subtract   ", AdditionSubtractionWorkload< boostmp_128_15_uint32_t >, NR_OPS / 2);
	PerformanceRunner("boostmp<256,15,uint32_t   add/subtract   ", AdditionSubtractionWorkload< boostmp_256_15_uint32_t >, NR_OPS / 4);
	PerformanceRunner("boostmp<512,15,uint32_t>  add/subtract   ", AdditionSubtractionWorkload< boostmp_512_15_uint32_t >, NR_OPS / 8);
	PerformanceRunner("boostmp<1024,15,uint32_t> add/subtract   ", AdditionSubtractionWorkload< boostmp_1024_15_uint32_t >, NR_OPS / 16);

	NR_OPS = 1024 * 32;
//	PerformanceRunner("boostmp<8,2,uint16_t>     division       ", DivisionWorkload< boostmp_8_2_uint8_t >, NR_OPS);
//	PerformanceRunner("boostmp<16,5,uint16_t>    division       ", DivisionWorkload< boostmp_16_5_uint16_t >, NR_OPS);
	PerformanceRunner("boostmp<32,8,uint32_t>    division       ", DivisionWorkload< boostmp_32_8_uint32_t >, NR_OPS);
	PerformanceRunner("boostmp<64,11,uint64_t>   division       ", DivisionWorkload< boostmp_64_11_uint64_t >, NR_OPS);
	PerformanceRunner("boostmp<128,15,uint32_t>  division       ", DivisionWorkload< boostmp_128_15_uint32_t >, NR_OPS / 2);
	PerformanceRunner("boostmp<256,15,uint32_t   division       ", DivisionWorkload< boostmp_256_15_uint32_t >, NR_OPS / 4);
	PerformanceRunner("boostmp<512,15,uint32_t>  division       ", DivisionWorkload< boostmp_512_15_uint32_t >, NR_OPS / 8);
	PerformanceRunner("boostmp<1024,15,uint32_t> division       ", DivisionWorkload< boostmp_1024_15_uint32_t >, NR_OPS / 16);

	// multiplication is the slowest operator

	NR_OPS = 1024 * 32;
//	PerformanceRunner("boostmp<8,2,uint16_t>     multiplication ", MultiplicationWorkload< boostmp_8_2_uint8_t >, NR_OPS);
	PerformanceRunner("boostmp<16,5,uint16_t>    multiplication ", MultiplicationWorkload< boostmp_16_5_uint16_t >, NR_OPS);
	PerformanceRunner("boostmp<32,8,uint32_t>    multiplication ", MultiplicationWorkload< boostmp_32_8_uint32_t >, NR_OPS);
	PerformanceRunner("boostmp<64,11,uint64_t>   multiplication ", MultiplicationWorkload< boostmp_64_11_uint64_t >, NR_OPS);
	PerformanceRunner("boostmp<128,15,uint32_t>  multiplication ", MultiplicationWorkload< boostmp_128_15_uint32_t >, NR_OPS / 2);
	PerformanceRunner("boostmp<256,15,uint32_t   multiplication ", MultiplicationWorkload< boostmp_256_15_uint32_t >, NR_OPS / 4);
	PerformanceRunner("boostmp<512,15,uint32_t>  multiplication ", MultiplicationWorkload< boostmp_512_15_uint32_t >, NR_OPS / 8);
	PerformanceRunner("boostmp<1024,15,uint32_t> multiplication ", MultiplicationWorkload< boostmp_1024_15_uint32_t >, NR_OPS / 16);

}

// conditional compilation
#define MANUAL_TESTING 0
#define STRESS_TESTING 0

int main()
try {
	using namespace sw::universal;

	std::string tag = "arithmetic performance benchmarking";

#if MANUAL_TESTING

	using Scalar = cfloat<128, 15, uint32_t>;
	Scalar a{ 1.0f }, b{ 1.0625f };;
	std::cout << a << " : " << b << std::endl;

	PerformanceRunner("cfloat<128,15,uint64_t>  multiplication ", MultiplicationWorkload< sw::universal::cfloat<128, 15, uint32_t> >, 1024);

	return 0;

	size_t NR_OPS = 10000000;
	PerformanceRunner("cfloat<16,5,uint16_t>    copy           ", CopyWorkload< sw::universal::cfloat<16, 5, uint16_t> >, NR_OPS);
	PerformanceRunner("cfloat<16,5,uint32_t>    copy           ", CopyWorkload< sw::universal::cfloat<16, 5, uint32_t> >, NR_OPS);


	std::cout << "done" << std::endl;

	return EXIT_SUCCESS;
#else
	std::cout << tag << std::endl;

	int nrOfFailedTestCases = 0;
	   
	TestCopyPerformance();
	TestDecodePerformance();
	TestArithmeticOperatorPerformance();
	TestBoostMPArithmeticOperatorPerformance();

	return (nrOfFailedTestCases > 0 ? EXIT_FAILURE : EXIT_SUCCESS);

#endif // MANUAL_TESTING
}
catch (char const* msg) {
	std::cerr << "Caught exception: " << msg << '\n';
	return EXIT_FAILURE;
}
catch (const std::runtime_error& err) {
	std::cerr << "Uncaught runtime exception: " << err.what() << std::endl;
	return EXIT_FAILURE;
}
catch (...) {
	std::cerr << "Caught unknown exception" << '\n';
	return EXIT_FAILURE;
}

/*
ETLO
Date run : 2/23/2020
Processor: Intel Core i7-7500 CPU @ 2.70GHz, 2 cores, 4 threads, 15W mobile processor
Memory   : 16GB
System   : 64-bit Windows 10 Pro, Version 1803, x64-based processor, OS build 17134.165

*/

/*
* ETLO
Date run : 1/22/2020
Processor: Intel Core i7-10700 CPU @ 2.90GHz, 8 cores, 16 threads, 45W desktop processor
Memory   : 16GB
System   : 64-bit Windows 10 Pro, Version 21H1, x64-based processor, OS build 19043.1466

arithmetic performance benchmarking
copy performance comparisoni
very large representations
cfloat<80,11,uint64_t>   copy              10000000 per           1e-07sec -> 100 Tops/sec
cfloat<96,11,uint64_t>   copy              10000000 per       0.0053574sec ->   1 Gops/sec
cfloat<128,11,uint64_t>  copy              10000000 per           1e-07sec -> 100 Tops/sec
cfloat<256,11,uint64_t>  copy              10000000 per           1e-07sec -> 100 Tops/sec
cfloat<512,11,uint64_t>  copy              10000000 per           1e-07sec -> 100 Tops/sec
classic floating-point cfloat decode operator performance
very large representations
cfloat<64,11,uint64_t>   decode              100000 per       0.0001097sec -> 911 Mops/sec
cfloat<80,11,uint64_t>   decode              100000 per        0.008986sec ->  11 Mops/sec
cfloat<96,11,uint64_t>   decode              100000 per       0.0106466sec ->   9 Mops/sec
cfloat<128,15,uint64_t>  decode              100000 per       0.0190613sec ->   5 Mops/sec
cfloat<256,15,uint64_t>  decode              100000 per       0.0302473sec ->   3 Mops/sec
cfloat<512,15,uint64_t>  decode              100000 per       0.0632031sec ->   1 Mops/sec
classic floating-point cfloat arithmetic operator performance
cfloat<8,2,uint8_t>      add/subtract       1000000 per       0.0035619sec -> 280 Mops/sec
cfloat<16,5,uint16_t>    add/subtract       1000000 per       0.0327843sec ->  30 Mops/sec
cfloat<32,8,uint32_t>    add/subtract       1000000 per       0.0329305sec ->  30 Mops/sec
cfloat<64,11,uint64_t>   add/subtract       1000000 per       0.0467038sec ->  21 Mops/sec
cfloat<128,11,uint32_t>  add/subtract        500000 per       0.0253118sec ->  19 Mops/sec
cfloat<128,15,uint32_t>  add/subtract        500000 per       0.0255983sec ->  19 Mops/sec
cfloat<256,15,uint32_t   add/subtract        250000 per       0.0211664sec ->  11 Mops/sec
cfloat<512,15,uint32_t>  add/subtract        125000 per       0.0158287sec ->   7 Mops/sec
cfloat<1024,15,uint32_t> add/subtract         62500 per       0.0129324sec ->   4 Mops/sec
cfloat<8,2,uint16_t>     division             32768 per       0.0041756sec ->   7 Mops/sec
cfloat<16,5,uint16_t>    division             32768 per       0.0091851sec ->   3 Mops/sec
cfloat<32,8,uint32_t>    division             32768 per       0.0249642sec ->   1 Mops/sec
cfloat<64,11,uint64_t>   division             32768 per       0.0358724sec -> 913 Kops/sec
cfloat<128,15,uint32_t>  division             16384 per        0.219809sec ->  74 Kops/sec
cfloat<256,15,uint32_t   division              8192 per        0.561821sec ->  14 Kops/sec
cfloat<512,15,uint32_t>  division              4096 per         1.35456sec ->   3 Kops/sec
cfloat<1024,15,uint32_t> division              2048 per         2.26855sec -> 902  ops/sec
cfloat<8,2,uint16_t>     multiplication       32768 per       0.0002174sec -> 150 Mops/sec
cfloat<16,5,uint16_t>    multiplication       32768 per       0.0015268sec ->  21 Mops/sec
cfloat<32,8,uint32_t>    multiplication       32768 per       0.0023463sec ->  13 Mops/sec
cfloat<64,11,uint64_t>   multiplication       32768 per       0.0002018sec -> 162 Mops/sec
cfloat<128,15,uint32_t>  multiplication       16384 per       0.0193484sec -> 846 Kops/sec
cfloat<256,15,uint32_t   multiplication        8192 per        0.037689sec -> 217 Kops/sec
cfloat<512,15,uint32_t>  multiplication        4096 per       0.0905458sec ->  45 Kops/sec
cfloat<1024,15,uint32_t> multiplication        2048 per        0.196303sec ->  10 Kops/sec
boost::multiprecision arithmetic operator performance
dummy case to fool the optimizer
boostmp<8,2,uint8_t>      add/subtract       1000000 per       0.0042658sec -> 234 Mops/sec
dummy case to fool the optimizer
boostmp<16,5,uint16_t>    add/subtract       1000000 per       0.0044779sec -> 223 Mops/sec
boostmp<32,8,uint32_t>    add/subtract       1000000 per       0.0344032sec ->  29 Mops/sec
boostmp<64,11,uint64_t>   add/subtract       1000000 per       0.0272771sec ->  36 Mops/sec
boostmp<128,11,uint32_t>  add/subtract        500000 per       0.0200895sec ->  24 Mops/sec
boostmp<128,15,uint32_t>  add/subtract        500000 per       0.0196782sec ->  25 Mops/sec
boostmp<256,15,uint32_t   add/subtract        250000 per       0.0115137sec ->  21 Mops/sec
boostmp<512,15,uint32_t>  add/subtract        125000 per       0.0081723sec ->  15 Mops/sec
boostmp<1024,15,uint32_t> add/subtract         62500 per       0.0069998sec ->   8 Mops/sec
boostmp<32,8,uint32_t>    division             32768 per       0.0014849sec ->  22 Mops/sec
boostmp<64,11,uint64_t>   division             32768 per       0.0040553sec ->   8 Mops/sec
boostmp<128,15,uint32_t>  division             16384 per       0.0044777sec ->   3 Mops/sec
boostmp<256,15,uint32_t   division              8192 per       0.0045808sec ->   1 Mops/sec
boostmp<512,15,uint32_t>  division              4096 per       0.0060263sec -> 679 Kops/sec
boostmp<1024,15,uint32_t> division              2048 per       0.0093697sec -> 218 Kops/sec
boostmp<16,5,uint16_t>    multiplication       32768 per       0.0001394sec -> 235 Mops/sec
boostmp<32,8,uint32_t>    multiplication       32768 per       0.0010901sec ->  30 Mops/sec
boostmp<64,11,uint64_t>   multiplication       32768 per       0.0006485sec ->  50 Mops/sec
boostmp<128,15,uint32_t>  multiplication       16384 per        0.000786sec ->  20 Mops/sec
boostmp<256,15,uint32_t   multiplication        8192 per       0.0007409sec ->  11 Mops/sec
boostmp<512,15,uint32_t>  multiplication        4096 per       0.0010394sec ->   3 Mops/sec
boostmp<1024,15,uint32_t> multiplication        2048 per       0.0017169sec ->   1 Mops/sec
*/
