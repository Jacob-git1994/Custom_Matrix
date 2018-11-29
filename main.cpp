#include "matrix.h"
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include "gauss.h"
//#include <omp.h>

/*Testing the capabilities of my matrix class with explict/implict heat equation*/
int main()
{
#if 1
	/*Set max threads*/
	//omp_set_num_threads(4);
	/*Initalizations*/
	FILE *f = fopen("data.txt","w");
	if(!f) {std::cout << "Error Opening File\n"; exit(1);} else {}
	const unsigned int size = 500;
	const unsigned int t_end = 1000000;
	nicole::matrix<double> grid = nicole::matrix<double>::linspace(0,1,size);
	nicole::matrix<double> inital = nicole::matrix<double>::zeros(size,1);
	nicole::matrix<double> time = nicole::matrix<double>::linspace(0,.5,t_end);
	
	/*Deletas*/
	const double dx = grid(2) - grid(1);
	const double dt = time(2) - time(1);
	
	/*Initalizations of the matrix we will use*/
	nicole::matrix<double> A = nicole::matrix<double>::zeros(size);
	nicole::matrix<double> b = nicole::matrix<double>::zeros(size,1);
	
	/*set bounary conditions*/
	const double l_bc = (double)1;
	const double r_bc = (double)1;
	
	/*Initalize matrix boundary data*/
	A(0,0) = (double)-2; b(0) = l_bc;
	A(0,1) = (double)1;
	A(size-1,size-1) = (double)-2; b(b.numel()-1) = r_bc;
	A(size-1,size-2) = (double)1;
	
#pragma omp parallel for simd 
	for(int i = 1; i < size-1; i++)
	{
		A(i,i-1) = (double)1;
		A(i,i) = (double)-2;
		A(i,i+1) = (double)1;
	}
	
	/*Redefine Matrix and function vector*/
	A *= (dt/std::pow(dx,2));
	b *= (dt/std::pow(dx,2));
	
	/*Start solving*/
	const unsigned int samples = 10000;
	unsigned int counter = samples;
	for(int i = 0; i <= time.numel(); i++)
	{
		if(counter == samples)
		{
			std::cout << (double)100*((double)i)/(double)time.numel() << "% done\n";
			for(int j = 0; j < size; j++)
			{
				fprintf(f,"%.14f\t%.14f\n",grid(j),inital(j));
			}
			counter = 0;
		}
		else
		{
			counter++;
		}
		inital += mat_mull_vec_band(A,inital,3) + b;
		//inital += A*inital + b; //Implict  = u^n+1  + A*u^n+1 (I + A)u^n+1 = u^n - b
	}
	fclose(f);
#else
	/*Set max threads*/
	//omp_set_num_threads(4);
	/*Initalizations*/
	FILE *f = fopen("data.txt","w");
	if(!f) {std::cout << "Error Opening File\n"; exit(1);} else {}
	const unsigned int size = 100;
	const unsigned int t_end = 100;
	nicole::matrix<double> grid = nicole::matrix<double>::linspace(0,1,size);
	nicole::matrix<double> inital = nicole::matrix<double>::zeros(size,1);
	nicole::matrix<double> time = nicole::matrix<double>::linspace(0,.5,t_end);
	
	/*Deletas*/
	const double dx = grid(2) - grid(1);
	const double dt = time(2) - time(1);
	
	/*Initalizations of the matrix we will use*/
	nicole::matrix<double> A = nicole::matrix<double>::zeros(size);
	nicole::matrix<double> b = nicole::matrix<double>::zeros(size,1);
	
	/*set bounary conditions*/
	const double l_bc = (double)1;
	const double r_bc = (double)0;
	
	/*Initalize matrix boundary data*/
	A(0,0) = (double)-2; b(0) = l_bc;
	A(0,1) = (double)1;
	A(size-1,size-1) = (double)-2; b(nicole::location::end) = r_bc;
	A(size-1,size-2) = (double)1;
	
#pragma omp parallel for simd 
	for(int i = 1; i < size-1; i++)
	{
		A(i,i-1) = (double)1;
		A(i,i) = (double)-2;
		A(i,i+1) = (double)1;
	}
	
	/*Redefine Matrix and function vector*/
	A *= (dt/std::pow(dx,2));
	b *= (dt/std::pow(dx,2));
	
	A -= nicole::matrix<double>::eye(A.rows());
	
	/*Start solving*/
	const unsigned int samples = 5;
	unsigned int counter = samples;
	for(int i = 0; i <= time.numel(); i++)
	{
		if(counter == samples)
		{
			std::cout << (double)100*((double)i)/(double)time.numel() << "% done\n";
			for(int j = 0; j < size; j++)
			{
				fprintf(f,"%.14f\t%.14f\n",grid(j),inital(j));
			}
			counter = 0;
		}
		else
		{
			counter++;
		}
		//nicole::matrix<double> temp = -(b + inital);
		//inital = nicole::gauss<double>::solve(A,-(b + inital));
		inital = nicole::jacobi<double>::solve(A,-(b + inital),.0001,10000);
	}
	//std::cout << "100% done!\n" << nicole::gauss<double>::det(A) << "\n";
	fclose(f);
#endif
}