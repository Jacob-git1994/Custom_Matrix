#include <stdlib.h>
#include "matrix.h"
#include "gauss.h"
#include <iostream>
//#include <omp.h>

using namespace nicole;
int main()
{
	//omp_set_num_threads(8);
	srand(time(NULL));
	const unsigned int N = 10;
	matrix<double> A = zeros<double>(N);
	
	A(0,0) = -(double)2;
	A(0,1) = (double)1;
	A(N-1,N-2) = (double)1;
	A(location::end,location::end) = (double)-2;
	
	for(int i = 1; i <N-1; i++)
	{
		A(i,i-1) = (double)1;
		A(i,i) = (double)-2;
		A(i,i+1) = (double)1;
	}
	
	std::cout << A << gauss<double>::det(A) << "\n";
	
}