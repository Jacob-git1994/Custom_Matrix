#include "matrix.h"
#include <iostream>
#include <omp.h>

using namespace nicole;
int main()
{
	const unsigned int N = 10000;
	mat A = dzeros(N);
	
	A(0,0) = -2.;
	A(0,1) = 1.;
	A(N-1,N-1) = -2.;
	A(N-1,N-2) = 1.;
	
#pragma omp parallel for simd
	for(int i = 1; i < N-1; i++)
	{
		A(i,i-1) = 1.;
		A(i,i) = -2.;
		A(i,i+1) = 1;
	}
	
	mat x = dlinspace(1,N,N).t();
	
	double start = omp_get_wtime();
	for(int i = 0; i < 10000; i++)
	{
		//A*x;
		mat_mull_vec_band(A,x,3);
	}
	double end = omp_get_wtime();	
	std::cout << end-start << "\n";
}