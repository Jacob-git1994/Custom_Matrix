#include <iostream>
#include "matrix.h"
#include "ode.h"
#include <cmath>

/*For the matrix class*/
using namespace nicole;

/*Grid and time constants*/
#define dt .0000001
#define dx .0000001

/*Constants to the fluid problem*/
#define h0 .01
#define eps h0/10.
#define mu 1.
#define sigma 1.
#define contact_angle 90.*pi/180.
#define n 3.
#define m 2.
#define h_star h0/100.

/*Class to store all the physical functions of the fluid*/
class fluid_functions
{
public:
	static double surface_tension(double,double,double,double,double);
	static double PIP(double);
	static double van_der_walls(double,double,double);
};

double fluid_functions::surface_tension(double hm2,double hm1,double h,double hp1,double hp2)
{
	return -sigma*(std::pow(hp1 + h,3)*(hp2 - 3.*hp1 + 3.*h - hm1) - std::pow(h + hm1,3)*(hp1 - 3.*h + 3.*hm1 - hm2)) \
		/ (mu*24.*std::pow(dx,4));
}

/*PI prime seperated to use later*/
double fluid_functions::PIP(double h)
{
	return (sigma*(std::cos(contact_angle) - 1.)*(m - 1.)*(n - 1.)*((h_star*m*std::pow((h_star/h),(m - 1.)))/std::pow(h,2.) - (h_star*n*std::pow((h_star/h),(n - 1.)))/std::pow(h,2.)))/(h_star*(m - n));
}

double fluid_functions::van_der_walls(double hm1,double h,double hp1)
{
	return -(std::pow(h + hp1,3)*fluid_functions::PIP((h + hp1)/2.)*(hp1 - h) - std::pow(h + hm1,3)*fluid_functions::PIP((h + hm1)/2.)*(h - hm1))/(24.*mu*std::pow(dx,2.));
}

/*Generate our function vector*/
static matrix<double> function_vector(const matrix<double> &v)
{
	/*Save the solution*/
	matrix<double> sol(v.rows(),1);
	
	/*Add in the boundary conditions*/
	/*Using Numan BC*/
	/*Left Boundary condition*/
	sol(0) = fluid_functions::surface_tension(v(1),v(0),v(0),v(1),v(2)) + fluid_functions::van_der_walls(v(0),v(0),v(1));
	sol(1) = fluid_functions::surface_tension(v(0),v(0),v(1),v(2),v(3)) + fluid_functions::van_der_walls(v(0),v(1),v(2));
	/*Right Boundary condition*/
	sol(v.rows()-1) = fluid_functions::surface_tension(v(v.rows()-3),v(v.rows()-2),v(v.rows()-1),v(v.rows()-1),v(v.rows()-2)) + fluid_functions::van_der_walls(v(v.rows()-2),v(v.rows()-1),v(v.rows()-1));
	sol(v.rows()-2) = fluid_functions::surface_tension(v(v.rows()-4),v(v.rows()-3),v(v.rows()-2),v(v.rows()-1),v(v.rows()-1)) + fluid_functions::van_der_walls(v(v.rows()-3),v(v.rows()-2),v(v.rows()-1));
	
	/*Interal grid*/
#pragma omp parallel for simd
	for(int i = 2; i < v.rows()-2; i++)
	{
		sol(i) = fluid_functions::surface_tension(v(i-2),v(i-1),v(i),v(i+1),v(i+2)) + fluid_functions::van_der_walls(v(i-1),v(i),v(i+1));
	}
	return sol;
}

int main()
{
	/*Set up our grid and spacing*/
	const matrix<double> grid = matrix<double>::delspace(0,.001,dx).t();
	const matrix<double> inital = h0 + eps*cos(2.*(double)pi*grid/grid(grid.numel()-1));
	
	/*solve*/
	ode<double>::rk4(inital,&function_vector,.01,dt,100);
	return 0;
}


