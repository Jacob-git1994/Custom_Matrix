#include <iostream>
#include "ode.h"
#include "matrix.h"
#include <cmath>
using namespace nicole;

/*Define constants*/
#define dt .0001
#define dx .05

/*Define boundary conditions*/
#define left_bc 1
#define right_bc 1

static double heat_equation(double hm1,double h,double hp1)
{
	return (hm1 - (double)2*h + hp1)/std::pow(dx,2);
}

static mat F(const mat &v)
{
	if(!v.is_col_vector()) {std::cout << "Not col vector\n"; exit(1);} else {}
	mat sol = dzeros(v.numel(),1);
	
	sol(0) = heat_equation(left_bc,v(0),v(1));
	sol(v.numel()) = heat_equation(v(v.numel()-2),v(v.numel()-1),right_bc);
	
#pragma omp parallel for simd
	for(int i = 1; i < v.numel()-1; i++)
	{
		sol(i) = heat_equation(v(i-1),v(i),v(i+1));
	}
	return sol;
}

int main()
{
	const mat grid = ddelspace(0,100,dx);
	const mat inital = sin(grid).t();
	ode<double>::rk4(inital,&F,200,dt,100);
}