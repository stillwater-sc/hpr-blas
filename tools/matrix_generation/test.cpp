// test structure to collaborate with Simunova
//
// Copyright (C) 2017-2018 Stillwater Supercomputing, Inc.
//
// Filename: eigenvalue_example.cpp (part of MTL4)

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

int main() {
	using namespace std;
	using Matrix = mtl::mat::dense2D< double >;
	Matrix M1(3, 3), M2(3, 3), M3(3, 3), M4(3, 3);

	M1 = 2, 0, 0,
		1, 1, 0,
		0, 1, 3; //EWs: 1,2,3     

	mtl::mat::eigenvalue_solver<Matrix> E1(M1);
	E1.setMaxIteration(10);
	E1.calc();
	cout << "M1(setting the number of iteraions): "
		<< E1.get_eigenvalues() << "\n";
}
