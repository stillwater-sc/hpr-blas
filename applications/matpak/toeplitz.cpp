#include <iostream>

#include <hprblas>
#include <matpak/toeplitz.hpp>


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
		using Vector = mtl::vec::dense_vector< Scalar >;
		Vector c{1, 2, 3, 4};
		Vector r{7, 5, 6, 7, 8};
		auto A = toeplitz(c,r);
		std::cout <<  A << std::endl;
 	}


	return 0;
}
