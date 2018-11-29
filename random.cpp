#include "random.hpp"
#include "matrix.h"
#include <iostream>
#include <iomanip>

int main()
{
	nicole::random_generator<double> gen(34012);
	nicole::matrix<double> X = gen.uniform((unsigned int)10,(unsigned int)10);
	std::cout << X;
	return 0;
}