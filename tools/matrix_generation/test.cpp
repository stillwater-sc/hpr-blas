// test structure to collaborate with Simunova
//
// Copyright (C) 2017-2018 Stillwater Supercomputing, Inc.
//
// Filename: eigenvalue_example.cpp (part of MTL4)

#include <hprblas>

using Matrix = mtl::mat::dense2D< double >;

int tutorial_example() {
	using namespace std;

	Matrix M1(3,3), M2(3,3), M3(3,3), M4(3,3);

	M1 = 2,0,0,
	     1,1,0,
	     0,1,3; //EWs: 1,2,3     

	mtl::mat::eigenvalue_solver<Matrix> E1(M1);
	E1.setMaxIteration(10);
	E1.calc();
	cout << "M1(setting the number of iterations): " 
	     << E1.get_eigenvalues() << "\n";

	M2 = 1,0,0,
	     0,1,5,
	     0,-2,3; //EWs: 1,2+3i,2-3i

	mtl::mat::eigenvalue_solver<Matrix> E2(M2);
	E2.setTolerance(1.0e-10);
	E2.calc();
	cout << "M2(providing tolerance): " 
	     << E2.get_eigenvalues() << "\n"; 

	M3 = -261, 209,  -49,
	     -530, 422,  -98,
	     -800, 631, -144; //EWs: 3,4,10  

	mtl::mat::eigenvalue_solver<Matrix> E3(M3);
	E3.setMaxIteration(10);
	E3.setTolerance(1.0e-10);
	E3.calc();
	cout << "M3(providing both): " 
	     << E3.get_eigenvalues() << "\n";

	M4 = 1,-3,3,
	     3,-5,3,
	     6,-6,4; //EWs: -2,-2,4

	mtl::mat::eigenvalue_solver<Matrix> E4(M4);
	E4.calc();   
	cout << "M4(with defaults): " 
	     << E4.get_eigenvalues() << "\n";    

	// Creating the solver implicitly
	cout << "M4(with defaults): " << eigenvalues(M4) << "\n";
}

int main() {
    using namespace std;
    using namespace sw::hprblas;

    tutorial_example();

    constexpr size_t N = 10;
    Matrix A(N,N);

    uniform_rand_diagonally_dominant(A);
    cout << eigenvalues(A) << endl;
}
