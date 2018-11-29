#include "pair.h"
#include "matrix.h"
#include "gauss.h"
#include <iostream>

using namespace nicole;

static unsigned int f(unsigned int x)
{
	return x*x;
}
int main()
{
	srand(time(NULL));
	matrix<unsigned int> x = matrix<unsigned int>::linspace(0,10,11);
	for(int i = 0; i < 2; i++)
	{
		x = row_concat(x,x);
	}
	std::cout << x.for_row_range(0,5,10,&f);
}