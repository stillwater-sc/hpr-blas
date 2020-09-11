
#include <iostream>
#include <hprblas>
#include <matpak/rowsto.hpp>


// Main function
int main()
{
    using namespace std;
	using namespace sw::unum;
	using namespace sw::hprblas;

    int m = 3;
    int n = 3;

    using Scalar = double;
    using Matrix = mtl::mat::dense2D< Scalar >;
	Matrix A = sw::hprblas::matpak::rowsto< Matrix >(m,n);   //
    cout <<  A << endl;
    // Matrix rowsto<Matrix>(m, n)

    return 0;
}
