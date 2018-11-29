#include "matrix.h"
#include "gauss.h"
#include <iostream>
#include <stdlib.h>

using namespace nicole;
int main()
{
	srand(time(NULL));
	
	const unsigned int N = 15;
	mat A = dzeros(N);
	mat b = drandu(N,1);
	A(0,0) = -2.; A(0,1) = 1.; A(0,2) = 1.;
	A(1,0) = 1.; A(1,1) = 2.; A(1,2) = 1.; A(1,3) = 1.;
	A(N-2,N-1) = 1.; A(N-2,N-2) = -2.; A(N-2,N-3) = 1.; A(N-2,N-4) = 1.;
	A(N-1,N-1) = -2.; A(N-1,N-2) = 1.; A(N-1,N-3) = 1.;
	
	for(int i = 2; i < N-2; i++)
	{
		A(i,i-2) = 1.;
		A(i,i-1) = 1.;
		A(i,i) = -2.;
		A(i,i+1) = 1.;
		A(i,i+2) = 1.;
	}
	
	std::cout << A << b << gauss<double>::solve_banded(A,b,5);
}