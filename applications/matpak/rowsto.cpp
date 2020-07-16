#include <iostream>

#include <hprblas>
#include <matpak/rowsto.hpp>

//#include <boost/numeric/mtl/mtl.hpp>

int main ()
{
	using namespace std;
	using namespace mtl;
	using namespace sw::unum;
	using namespace sw::hprblas;
	cout << setprecision(20);
	{
		using Scalar = double;
		using Matrix = mtl::mat::dense2D< Scalar >;
		Matrix A = sw::hprblas::matpak::rowsto< Matrix >(5,5);   //
		std::cout <<  A << std::endl;
		for (size_t i=0;i<5;++i){
			cout << sum(A[i][iall])<<endl;
		}
	}

	{
		constexpr size_t nbits = 32;
		constexpr size_t es = 2;

		using Scalar = posit<nbits, es>;
		using Matrix = mtl::mat::dense2D< Scalar >;
		Matrix A = sw::hprblas::matpak::rowsto< Matrix >(5,5);   //
		std::cout <<  A << std::endl;
		for (size_t i=0;i<5;++i){
			cout << sum(A[i][iall])<<endl;
		}
	}

	return 0;
}
