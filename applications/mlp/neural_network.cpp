// neural_network.cpp example program using posit arithmetic for a neural network simulation
//
// Copyright (C) 2017-2018 Stillwater Supercomputing, Inc.
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.
#include "common.hpp"
// include the posit number system
#include <universal/number/posit/posit>
#define MTL_WITH_INITLIST 1
#include <boost/numeric/mtl/mtl.hpp>

#undef NOW
#ifdef NOW
template<size_t nbits, size_t es>
sw::universal::posit<nbits, es> Sigmoid_(sw::universal::posit<nbits, es>& x, bool derivative = false) {
	if (derivative) return x*(1-x);
	return (1);
}


template<typename Scalar>
Scalar sigmoid(Scalar x) {
	return (1.0f / (1.0f + exp(-x)));
}

template<typename Vector, typename Scalar>
Scalar f_theta(Scalar x, Scalar b, const Vector& V, const Vector& W) {
	double result = b;
	for (int i = 0; i < N; ++i) {
		result += V[i] * sigmoid(c[i] + W[i] * x);
	}
	return result;
}

template<typename Vector, typename Scalar>
void train(Scalar x, Scalar y, Vector& W, Vector& V, Vector& c) {
	for (int i = 0; i < size(W); ++i) {
		W[i] = W[i] - epsilon * 2 * (f_theta(x) - y) * V[i] * x *
			(1 - sigmoid(c[i] + W[i] * x)) * sigmoid(c[i] + W[i] * x);
	}
	for (int i = 0; i < size(V); ++i) {
		V[i] = V[i] - epsilon * 2 * (f_theta(x) - y) * sigmoid(c[i] + W[i] * x);
	}
	b = b - epsilon * 2 * (f_theta(x) - y);
	for (int i = 0; i < size(c); ++i) {
		c[i] = c[i] - epsilon * 2 * (f_theta(x) - y) * V[i] *
			(1 - sigmoid(c[i] + W[i] * x)) * sigmoid(c[i] + W[i] * x);
	}
}

// Multi-Level Perceptron
// one input layer
// one hidden layer with NrOfNeurons. 
// A \(sigmoid\\) function as activation function
// one output layer
int main(int argc, char** argv)
try {
	using namespace std;
	using namespace mtl;
	using namespace sw::universal;

	const size_t nbits = 16;
	const size_t es = 1;
	const size_t vecSize = 32;

	using Scalar = float;
	using Matrix = mtl::dense2D<Scalar>;
	using Vector = mtl::dense_vector<Scalar>;
	using Pair   = std::pair<Scalar, Scalar>;
	using VoP    = mtl::dense_vector<Pair>;

#if defined(MTL_WITH_INITLIST) && defined(MTL_WITH_AUTO) && defined(MTL_WITH_RANGEDFOR)
	Matrix X = {
		{0,0,1},
		{0,1,1},
		{1,0,1},
		{1,1,1}
	};
#endif

	constexpr int NrOfTrainingSamples = 20;
	constexpr int NrOfNeurons = 5;
	constexpr float epsilon = 0.05f;
	constexpr int epoch = 50000;

	Vector c(NrOfNeurons), W(NrOfNeurons), V(NrOfNeurons);
	c = W = V = 0;
	Scalar b = 0;

	// fill initial state with random values
	srand((unsigned int)time(NULL));
	for (int i = 0; i < NrOfNeurons; i++) {
		W[i] = Scalar(2 * rand() / RAND_MAX - 1);
		V[i] = Scalar(2 * rand() / RAND_MAX - 1);
		c[i] = Scalar(2 * rand() / RAND_MAX - 1);
	}
	VoP trainingSet(NrOfTrainingSamples);

	for (int i = 0; i < NrOfTrainingSamples; i++) {
		trainingSet[i] = make_pair(i * 2 * m_pi / NrOfTrainingSamples, sin(i * m_2_pi / NrOfTrainingSamples));
	}

	// train for epoch number of cycles
	for (int k = 0; k < epoch; ++k) {
		for (int i = 0; i < NrOfTrainingSamples; i++) {
			auto x = trainingSet[i].first;
			auto y = trainingSet[i].second;
			train(x,y, W, V, c);
		}
		std::cout << k << "\r";
	}

	// Plot the results
	int nrVisualSamples = 1000;
	Vector x(nrVisualSamples);
	Vector y1(nrVisualSamples), y2(nrVisualSamples);

	auto scaling_factor = m_2pi / nrVisualSamples;
	for (int i = 0; i < nrVisualSamples; i++) {
		x[i] =  i * scaling_factor;
		y1[i] = sin(i * scaling_factor);
		y2[i] = f_theta(Scalar(i * scaling_factor), b, V, W);
	}

	FILE * gp = _popen("gnuplot", "w");
	fprintf(gp, "set terminal wxt size 600,400 \n");
	fprintf(gp, "set grid \n");
	fprintf(gp, "set title '%s' \n", "f(x) = sin (x)");
	fprintf(gp, "set style line 1 lt 3 pt 7 ps 0.1 lc rgb 'green' lw 1 \n");
	fprintf(gp, "set style line 2 lt 3 pt 7 ps 0.1 lc rgb 'red' lw 1 \n");
	fprintf(gp, "plot '-' w p ls 1, '-' w p ls 2 \n");

	// Exact f(x) = sin(x) -> Green Graph
	for (int k = 0; k < nrVisualSamples; ++k) {
		fprintf(gp, "%f %f \n", x[k], y1[k]);
	}
	fprintf(gp, "e\n");

	//Neural Network Approximate f(x) = sin(x) -> Red Graph
	for (int k = 0; k < nrVisualSamples; ++k) {
		fprintf(gp, "%f %f \n", x[k], y2[k]);
	}
	fprintf(gp, "e\n");

	fflush(gp);

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

#else
int main() {}
#endif
