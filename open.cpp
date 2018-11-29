#include "matrix.h"
#include "open_blas.hpp"
#include <stdio.h>

using namespace nicole;
int main()
{	
	//10x5 * 5x6 = 10x6 + 10x6
	const unsigned int N = 10000;
	const mat A = drandu(0,10,N);
	const mat x = drandu(0,10,N);
	
	std::cout << cblas<double>::cblas_gemm(A,x,dzeros(N,N),1.,0.);
	
}