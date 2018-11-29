#include <omp.h>
#include "matrix.h"
#include "gauss.h"
#include <iostream>
#include <stdlib.h>

using namespace nicole;
int main()
{
	srand(time(NULL));
	double start,end,begin;
	unsigned int size = 2000;
	
	for(int i = 1; i < 5; i++)
	{
		omp_set_num_threads(i);
		matrix<double> A = matrix<double>::randu(0,100,size);
		start = omp_get_wtime();
		matrix<double> B = A*A;
		end = omp_get_wtime();
		if(i == 1) {begin = end-start;} else {}
		std::cout << "Tp: " << end-start << " Sp: " << begin/(end-start) <<  std::endl;
	}
}