// Reflect Matrix across counter-identity

#include <iostream>
#include <hprblas>
#include <matpak/rowsto.hpp>
#include <matpak/fliptranspose.hpp>



int main ()
{
	using namespace std;
	using namespace mtl;
	using namespace sw::unum;
	using namespace sw::hprblas;
	using namespace sw::hprblas::matpak;
	cout << setprecision(8);
	{
		using Scalar = posit<16,1>;
		// using Scalar = double;

		using Matrix = mtl::mat::dense2D< Scalar >;
		Matrix A = rowsto< Matrix >(5,5);   //
		std::cout <<  A << std::endl;
		std::cout <<  fliptranspose(A) << std::endl;
 	}

	return 0;
}
