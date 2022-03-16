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

#undef MP_EXPERIMENT
#ifdef MP_EXPERIMENT
#define UNIVERSAL_CONSTEXPR constexpr 
#define UNIVERSAL_STATIC_ASSERT( ... ) static_assert(__VA_ARGS__, #__VA_ARGS__)
// UNIVERSAL_STATIC_CONSTANT workaround --------------------------------------- //
// On compilers which don't allow in-class initialization of static integral
// constant members, we must use enums as a workaround if we want the constants
// to be available at compile-time. This macro gives us a convenient way to
// declare such constants.

#  ifdef UNIVERSAL_NO_INCLASS_MEMBER_INITIALIZATION
#       define UNIVERSAL_STATIC_CONSTANT(type, assignment) enum { assignment }
#  else
#     define UNIVERSAL_STATIC_CONSTANT(type, assignment) static const type assignment
#  endif

namespace sw {
	namespace experimental {

		// rebind.hpp
//		namespace backends {
			namespace detail
			{
				template <class value_type, class my_allocator>
				struct rebind
				{
					typedef typename std::allocator_traits<my_allocator>::template rebind_alloc<value_type> type;
				};
			}
//		}

		// integral_c_tag.hpp
		struct integral_c_tag {
			static const int value = 0;
		};
		// integral constant.hpp

		namespace mpl
		{
			template <bool B> struct bool_;
			template <class I, I val> struct integral_c;
			struct integral_c_tag;
		}

		template <class T, T val>
		struct integral_constant
		{
			typedef integral_c_tag tag;
			typedef T value_type;
			typedef integral_constant<T, val> type;
			static const T value = val;

			operator const mpl::integral_c<T, val>& ()const
			{
				static const char data[sizeof(long)] = { 0 };
				static const void* pdata = data;
				return *(reinterpret_cast<const mpl::integral_c<T, val>*>(pdata));
			}
			UNIVERSAL_CONSTEXPR operator T()const { return val; }
		};

		template <class T, T val>
		T const integral_constant<T, val>::value;

		template <bool val>
		struct integral_constant<bool, val>
		{
			typedef integral_c_tag tag;
			typedef bool value_type;
			typedef integral_constant<bool, val> type;
			static const bool value = val;

			operator const mpl::bool_<val>& ()const
			{
				static const char data[sizeof(long)] = { 0 };
				static const void* pdata = data;
				return *(reinterpret_cast<const mpl::bool_<val>*>(pdata));
			}
			BOOST_CONSTEXPR operator bool()const { return val; }
		};

		template <bool val>
		bool const integral_constant<bool, val>::value;

		typedef integral_constant<bool, true> true_type;
		typedef integral_constant<bool, false> false_type;

		// is_void.hpp
		template <class T>
		struct is_void : public false_type {};

		template<> struct is_void<void> : public true_type {};
		template<> struct is_void<const void> : public true_type {};
		template<> struct is_void<const volatile void> : public true_type {};
		template<> struct is_void<volatile void> : public true_type {};


		using limb_type = unsigned long;

		enum integer_type {
			signed_magnitude = 0,
			unsigned_magnitude = 1
		};

		template<unsigned nbits>
		class number {
		public:

		};

		template<unsigned nbits, integer_type SignType, class Allocator>
		class universal_int_base
		{
			typedef typename detail::rebind<limb_type, Allocator>::type            allocator_type;
			typedef typename std::allocator_traits<allocator_type>::pointer        limb_pointer;
			typedef typename std::allocator_traits<allocator_type>::const_pointer  const_limb_pointer;
			UNIVERSAL_STATIC_ASSERT(!is_void<Allocator>::value);
		};

		template<unsigned nbits, integer_type SignType, class Allocator = void>
		class universal_int_backend 
			: public universal_int_base<nbits, SignType, Allocator>
		{
			typedef universal_int_backend<nbits, SignType, Allocator>  self_type;
		public:
		private:

			unsigned m_limbs;
//			bool     m_sign, m_internal;
		};

//		typedef number<universal_int_backend<> >   cpp_int;
}}
#endif // MP_EXPERIMENT

// measure performance of arithmetic operators
void TestUniversalIntegerOperatorPerformance() {
	using namespace sw::universal;
	std::cout << "\nUniversal Integer Arithmetic operator performance\n";

	uint64_t NR_OPS = 1000000;

	PerformanceRunner("integer<8>    add/subtract  ", AdditionSubtractionWorkload< sw::universal::integer<8,  uint8_t > >, NR_OPS);
	PerformanceRunner("integer<16>   add/subtract  ", AdditionSubtractionWorkload< sw::universal::integer<16, uint16_t> >, NR_OPS);
	PerformanceRunner("integer<32>   add/subtract  ", AdditionSubtractionWorkload< sw::universal::integer<32, uint32_t> >, NR_OPS);
	PerformanceRunner("integer<64>   add/subtract  ", AdditionSubtractionWorkload< sw::universal::integer<64, uint64_t> >, NR_OPS);
	PerformanceRunner("integer<128>  add/subtract  ", AdditionSubtractionWorkload< sw::universal::integer<128, uint32_t> >, NR_OPS / 2);
	PerformanceRunner("integer<256>  add/subtract  ", AdditionSubtractionWorkload< sw::universal::integer<256, uint32_t> >, NR_OPS / 4);
	PerformanceRunner("integer<512>  add/subtract  ", AdditionSubtractionWorkload< sw::universal::integer<512, uint32_t> >, NR_OPS / 8);
	PerformanceRunner("integer<1024> add/subtract  ", AdditionSubtractionWorkload< sw::universal::integer<1024, uint32_t> >, NR_OPS / 16);
	return;
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
	std::cout << "\nBoost Integer Arithmetic operator performance\n";


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

#undef SIZE_EXPERIMENT
#ifdef SIZE_EXPERIMENT
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

#undef MP_EXPERIMENT
#ifdef MP_EXPERIMENT
	sw::experimental::universal_int_backend<32, sw::experimental::signed_magnitude, void> a;
	cout << sizeof(a) << endl;
#endif // MP_EXPERIMENT

#ifdef CARRY_PROPAGATION_TEST
	// test to see how Boost handles carry propagation among limbs
	int1024_t a, b, c;
	a = 0xffff'ffff'ffff'ffff;
	b = (a << 64) + a;
	b = (b << 64) + a;
	b = (b << 64) + a;
	c = a + b;
	std::cout << a << " + " << b << " = " << c << '\n';
#endif

//	uint64_t NR_OPS = 100'000'000;
//  PerformanceRunner("integer<1024> add/subtract", AdditionSubtractionWorkload< sw::universal::integer<1024, uint32_t> >, NR_OPS);
//	PerformanceRunner("integer<128>  add/subtract  ", AdditionSubtractionWorkload< int128_t >, NR_OPS);
//	PerformanceRunner("integer<256>  add/subtract  ", AdditionSubtractionWorkload< int256_t >, NR_OPS);
//	PerformanceRunner("integer<512>  add/subtract  ", AdditionSubtractionWorkload< int512_t >, NR_OPS);
//	PerformanceRunner("integer<1024> add/subtract  ", AdditionSubtractionWorkload< int1024_t >, NR_OPS);

//	PerformanceRunner("cfloat<128>   multiplication", MultiplicationWorkload< sw::universal::cfloat<128,15,uint32_t, true, false, false> >, NR_OPS / 8);
//	PerformanceRunner("integer<128>  multiplication", MultiplicationWorkload< sw::universal::integer<128,uint32_t> >, NR_OPS / 8);
//	PerformanceRunner("integer<512>  multiplication", MultiplicationWorkload< sw::universal::integer<512,uint32_t> >, NR_OPS / 16);
//	PerformanceRunner("integer<1024> multiplication", MultiplicationWorkload< sw::universal::integer<1024,uint32_t> >, NR_OPS / 32);

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

/* performance comparison: 3/14/2022
Universal Integer Arithmetic operator performance
integer<8>    add/subtract      1000000 per       0.0007622sec ->   1 Gops/sec
integer<16>   add/subtract      1000000 per       0.0078274sec -> 127 Mops/sec
integer<32>   add/subtract      1000000 per       0.0096869sec -> 103 Mops/sec
integer<64>   add/subtract      1000000 per       0.0136763sec ->  73 Mops/sec
integer<128>  add/subtract       500000 per       0.0123696sec ->  40 Mops/sec
integer<256>  add/subtract       250000 per       0.0088511sec ->  28 Mops/sec
integer<512>  add/subtract       125000 per       0.0087157sec ->  14 Mops/sec
integer<1024> add/subtract        62500 per       0.0085516sec ->   7 Mops/sec
integer<8>    division            32768 per       0.0025333sec ->  12 Mops/sec
integer<16>   division            32768 per       0.0041453sec ->   7 Mops/sec
integer<32>   division            32768 per       0.0076186sec ->   4 Mops/sec
integer<64>   division            16384 per        0.006965sec ->   2 Mops/sec
integer<128>  division             8192 per       0.0062281sec ->   1 Mops/sec
integer<512>  division             4096 per        0.012138sec -> 337 Kops/sec
integer<1024> division             2048 per       0.0122465sec -> 167 Kops/sec
integer<8>    multiplication      32768 per       0.0055092sec ->   5 Mops/sec
integer<16>   multiplication      32768 per        0.017206sec ->   1 Mops/sec
integer<32>   multiplication      16384 per        0.035561sec -> 460 Kops/sec
integer<64>   multiplication       8192 per       0.0638777sec -> 128 Kops/sec
integer<128>  multiplication       4096 per        0.127786sec ->  32 Kops/sec
integer<512>  multiplication       2048 per         1.05581sec ->   1 Kops/sec
integer<1024> multiplication       1024 per         2.01186sec -> 508  ops/sec

Boost Integer Arithmetic operator performance
integer<128>  add/subtract       500000 per       0.0056359sec ->  88 Mops/sec
integer<256>  add/subtract       250000 per       0.0027532sec ->  90 Mops/sec
integer<512>  add/subtract       125000 per       0.0015111sec ->  82 Mops/sec
integer<1024> add/subtract        62500 per        0.000767sec ->  81 Mops/sec
integer<128>  division             8192 per       0.0001457sec ->  56 Mops/sec
integer<512>  division             4096 per        9.09e-05sec ->  45 Mops/sec
integer<1024> division             2048 per        3.66e-05sec ->  55 Mops/sec
integer<128>  multiplication       4096 per        4.98e-05sec ->  82 Mops/sec
integer<512>  multiplication       2048 per        2.96e-05sec ->  69 Mops/sec
integer<1024> multiplication       1024 per        2.73e-05sec ->  37 Mops/sec


3/16/2022: after introducing BlockType and pointer-based block visitation

Universal Integer Arithmetic operator performance
integer<8>    add/subtract      1000000 per       0.0005621sec ->   1 Gops/sec
integer<16>   add/subtract      1000000 per       0.0005484sec ->   1 Gops/sec
integer<32>   add/subtract      1000000 per       0.0004982sec ->   2 Gops/sec
integer<64>   add/subtract      1000000 per       0.0003297sec ->   3 Gops/sec
integer<128>  add/subtract       500000 per       0.0048675sec -> 102 Mops/sec
integer<256>  add/subtract       250000 per       0.0025151sec ->  99 Mops/sec
integer<512>  add/subtract       125000 per       0.0017968sec ->  69 Mops/sec
integer<1024> add/subtract        62500 per       0.0016404sec ->  38 Mops/sec

Boost Integer Arithmetic operator performance
integer<128>  add/subtract       500000 per       0.0052232sec ->  95 Mops/sec
integer<256>  add/subtract       250000 per       0.0026663sec ->  93 Mops/sec
integer<512>  add/subtract       125000 per       0.0013526sec ->  92 Mops/sec
integer<1024> add/subtract        62500 per       0.0007275sec ->  85 Mops/sec
integer<128>  division             8192 per       0.0001588sec ->  51 Mops/sec
integer<512>  division             4096 per        7.28e-05sec ->  56 Mops/sec
integer<1024> division             2048 per        3.94e-05sec ->  51 Mops/sec
integer<128>  multiplication       4096 per        4.92e-05sec ->  83 Mops/sec
integer<512>  multiplication       2048 per        2.68e-05sec ->  76 Mops/sec
integer<1024> multiplication       1024 per        1.51e-05sec ->  67 Mops/sec


3/16/2022: identified that boost performance comes from the fact that it is
dynamic in its limbs and the benchmark was feeding it single limb values.
That explains why when you go to larger integer sizes the performance didn't
drop as expected. Universal on the other hand will always do all the limb
work. So here is the boost performance when you give it full limbs:

Boost Integer Arithmetic operator performance
integer<128>  add/subtract      1000000 per        0.017574sec ->  56 Mops/sec
integer<256>  add/subtract      1000000 per       0.0224437sec ->  44 Mops/sec
integer<512>  add/subtract      1000000 per       0.0310884sec ->  32 Mops/sec
integer<1024> add/subtract      1000000 per       0.0590956sec ->  16 Mops/sec

Universal Integer Arithmetic operator performance
integer<128>  add/subtract       500000 per       0.0048675sec -> 102 Mops/sec
integer<256>  add/subtract       250000 per       0.0025151sec ->  99 Mops/sec
integer<512>  add/subtract       125000 per       0.0017968sec ->  69 Mops/sec
integer<1024> add/subtract        62500 per       0.0016404sec ->  38 Mops/sec

required this workload function and a mod to PerformanceRunner

	// special workload for adaptive integer types, like adaptint and boost::multiprecision
	// DO NOT USE as it requires manual mods to the performance runner
	template<typename Scalar>
	void AdditionSubtractionAdaptiveIntegerWorkload(size_t nbits) {
		size_t NR_OPS = 1000'000;
		Scalar a = 0x00FF'FFFF;
//		std::cout << a << '\n';
		a <<= (nbits - 32);
		Scalar b = a + 1;
//		std::cout << a << '\n';
//		std::cout << b << '\n';
		std::vector<Scalar> data = { a, -b };
		for (size_t i = 1; i < NR_OPS; ++i) {
			a = data[i % 2];
			b = b + a;
		}
		if (b == Scalar(0.0f)) {
			std::cout << "dummy case to fool the optimizer\n";
		}
	}

	PerformanceRunner("integer<128>  add/subtract  ", AdditionSubtractionAdaptiveIntegerWorkload< int128_t >, 128);
	PerformanceRunner("integer<256>  add/subtract  ", AdditionSubtractionAdaptiveIntegerWorkload< int256_t >, 256);
	PerformanceRunner("integer<512>  add/subtract  ", AdditionSubtractionAdaptiveIntegerWorkload< int512_t >, 512);
	PerformanceRunner("integer<1024> add/subtract  ", AdditionSubtractionAdaptiveIntegerWorkload< int1024_t >, 1024);

*/