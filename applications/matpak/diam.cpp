// Calculate the diameter of a matrix.
// Functions needed: max, norm, mtl::iall=:
#include <iostream>
#include <hprblas>
#include <matpak/ek.hpp>



int main ()
{
	using namespace std;
	using namespace mtl;
	using namespace sw::unum;
	using namespace sw::hprblas;
	using namespace sw::hprblas::matpak;
	cout << setprecision(5);
	{
		using Scalar = posit<16,1>;
		// using Scalar = double;
		auto v = ek<Scalar>(3,8);
		std::cout <<  v << std::endl;
 	}

	return 0;
}
