#include <iostream>

#include <hprblas>
#include <matpak/rowsto.hpp>
#include <matpak/rot90.hpp>


int main ()
{
	using namespace std;
	using namespace mtl;
	using namespace sw::unum;
	using namespace sw::hprblas;
	using namespace sw::hprblas::matpak;
	cout << setprecision(5);
	{
		using Scalar = double;
		using Matrix = mtl::mat::dense2D< Scalar >;
		Matrix A = rowsto< Matrix >(5,5);   //
		std::cout <<  A << std::endl;
		//auto B = rot90(A);
		//auto B = trans(A);
		//std::cout <<  B << std::endl;
	}

	{
		constexpr size_t nbits = 32;
		constexpr size_t es = 2;

		using Scalar = posit<nbits, es>;
		using Matrix = mtl::mat::dense2D< Scalar >;
		//Matrix A = rowsto< Matrix >(5,5);   //
		//std::cout <<  A << std::endl;
		//fliplr(A);
		//std::cout <<  A << std::endl;
	}

	return 0;
}
