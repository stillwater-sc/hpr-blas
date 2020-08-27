#include <iostream>

#include <hprblas>
#include <matpak/rowsto.hpp>
#include <matpak/isnormal.hpp>


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
		Matrix A_transpose(mtl::mat::trans(A));
		std::cout <<  A_transpose << std::endl;
		std::cout <<  A*A_transpose << std::endl;
		std::cout <<  A_transpose*A << std::endl;

		if(isnormal(A,0.00001)){
			std::cout <<  "Matrix A is normal" << std::endl;
		}else{
			std::cout <<  "Matrix A is NOT normal" << std::endl;
		}

	}

	{
		constexpr size_t nbits = 32;
		constexpr size_t es = 2;

		using Scalar = posit<nbits, es>;
		using Matrix = mtl::mat::dense2D< Scalar >;

	}

	return 0;
}
