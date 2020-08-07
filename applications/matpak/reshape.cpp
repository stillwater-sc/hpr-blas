#include <iostream>
#include <hprblas>
#include <matpak/rowsto.hpp>
#include <matpak/reshape.hpp>


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
		typedef mtl::mat::parameters<col_major, index::c_index, mtl::non_fixed::dimensions, false, unsigned long> column_matrix;
		using Matrix = mtl::dense2D<Scalar, column_matrix > ;
		Matrix A = rowsto< Matrix >(6,2);   //
		std::cout <<  A << std::endl;
		auto T = reshape(A,4,3);
		std::cout <<  T << std::endl;
	}
	return 0;
}
