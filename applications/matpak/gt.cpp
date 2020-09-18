#include <iostream>

#include <hprblas>
#include <matpak/rowsto.hpp>
#include <matpak/gt.hpp>

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
		Matrix B = rowsto< Matrix >(5,5);   //

		// B = A;

		std::cout <<  gt(A, B) << std::endl;
		auto C = A > B;
		std::cout <<  C << std::endl;
	}	

	
	return 0;
}
