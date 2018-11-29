#include "laplace.hpp"
#include <iostream>
#include <cmath>
using namespace nicole;

double forceing(double x,double y)
{
	return 0.;//std::sin(pi*x/10.);
}

double b(double x)
{
	return 1.;
}

double de(double x)
{
	return 0.;
}

int main()
{
	laplace<double> v(60,0,1);
	v.set_forcing(&forceing);
	v.set_boundary(&de,&de,&b,&b);
	v.set_boundary_type('d','d','d','d');
	v.initalize();
	v.lu_solve();
}
