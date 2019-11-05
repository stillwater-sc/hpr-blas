// finite_difference.cpp example program to demonstrate finite difference calculations
//
// Copyright (C) 2017-2019 Stillwater Supercomputing, Inc.
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.
#include "common.hpp"
#include <boost/multiprecision/cpp_int.hpp>
#define POSIT_THROW_ARITHMETIC_EXCEPTION 1
#include <universal/posit/posit>

/* concepts
   The complete definition of a template function or class should contain the list of required
   concepts as is done for functions in the Standard template Library; see http://www.sgi.com/tech/stl/

   */
template <typename F, typename T>
class forward_difference {
public:
	forward_difference(const F& f, const T& h) : f(f), h(h) {}
	T operator()(const T& x) const { return (f(x+h) - f(x)) / h; }
private:
	const F& f;
	T        h;
};

template<typename F, typename T>
class backward_difference {
public:
	backward_difference(const F& f, const T& h) : f(f), h(h) {}
	T operator()(const T& x) const { return (f(x) - f(x - h)) / h; }
private:
	const F& f;
	T        h;
};

template<typename F, typename T>
class central_difference {
public:
	central_difference(const F& f, const T& h) : f(f), h(h) {}
	T operator()(const T& x) const { return (f(x + T(0.5)*h) - f(x - T(0.5)*h)) / h; }
private:
	const F& f;
	T        h;
};

// Use recursion to define arbitrary n-th order derivatives

template<unsigned N, typename F, typename T>
class nth_derivative {
	using n_minus_1_derivative = nth_derivative<N - 1, F, T>;
public:
	nth_derivative(const F& f, const T& h) : h(h), fprev(f, h) {}
	T operator()(const T& x) const { return (fprev(x + h) - fprev(x)) / h; }
private:
	n_minus_1_derivative fprev;
	T                    h;
};

// termination condition
template<typename F, typename T>
class nth_derivative<1, F, T> : public forward_difference<F, T> {
	//using derivative<F, T>::derivative;
};

// trigonometry sin functor
template<typename T>
struct sin_f {
	T operator() (const T& x) const { return sin(x); }
};

// a functor with a scale parameter
template<typename T>
class scaled_sin_f {
public:
	scaled_sin_f(const T& scale) : scale(scale) {}
	T operator() (const T& x) const { return scale * sin(x); }
private:
	T scale;
};

// trigonometry cos functor
template<typename T>
struct cos_f {
	T operator() (const T& x) const { return cos(x); }
};

// trigonometry tan functor
template<typename T>
struct tan_f {
	T operator() (const T& x) const { return tan(x); }
};

template<typename Real>
void EnumerateDerivativeError() {
	using namespace std;
	using namespace sw::unum; // for m_pi_4

	cout << "sin at PI/4 : " << sin_f<Real>()(Real(m_pi_4)) << endl;
	cout << "cos at PI/4 : " << cos_f<Real>()(m_pi_4) << endl;
	cout << "tan at PI/4 : " << tan_f<Real>()(m_pi_4) << endl;

	using forward_sin_f = forward_difference<sin_f<Real>, Real>;
	using backward_sin_f = backward_difference<sin_f<Real>, Real>;
	using central_sin_f = central_difference<sin_f<Real>, Real>;

	constexpr int WIDTH = 16;

	Real h = Real(0.05);
	sin_f<Real> sin_o;
	cout << setw(WIDTH) << "step (h)" << setw(WIDTH) << "backward" << setw(WIDTH) << "error" << setw(WIDTH) << "forward" << setw(WIDTH) << "error" << setw(WIDTH) << "central" << setw(WIDTH) << "error\n";
	for (unsigned i = 0; i < 10; ++i) {
		// backward = left
		backward_sin_f ld_sin_o(sin_o, h);
		Real ld_sin = ld_sin_o(Real(m_pi_4));
		Real lerror = ld_sin - cos(m_pi_4);

		// forward = right
		forward_sin_f rd_sin_o(sin_o, h);
		Real rd_sin = rd_sin_o(Real(m_pi_4));
		Real rerror = rd_sin - cos(m_pi_4);

		// central = middle
		central_sin_f md_sin_o(sin_o, h);
		Real md_sin = md_sin_o(Real(m_pi_4));
		Real merror = md_sin - cos(m_pi_4);

		cout << setw(WIDTH) << h << setw(WIDTH) << ld_sin << setw(WIDTH) << lerror << setw(WIDTH) << rd_sin << setw(WIDTH) << rerror << setw(WIDTH) << md_sin << setw(WIDTH) << merror << endl;

		h *= 0.5;
	}
}

int main() 
try {
	using namespace std;
	using namespace boost::multiprecision;
	using namespace sw::unum;

	EnumerateDerivativeError<float>();

	EnumerateDerivativeError< posit<32, 2> >();
	
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
