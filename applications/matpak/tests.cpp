#include <iostream>
#include <hprblas>

using namespace std;

//using namespace sw::unum;
//using namespace hprblas;

void incre(int& x)
{
	x++;
}

int main ()
{
	int i1 = 2;
	int i2 = 5;
	float pi = 3.1415;
	double x = -1.5e6;
	char c1 = 'a';
	bool cmp = true;

	incre(i1);   // post successor, ++i1 increment first
	std::cout <<  i1 << std::endl;
	return 0;
}
