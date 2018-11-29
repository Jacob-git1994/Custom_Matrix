#include "matrix.h"
#include "conjgrad.hpp"
#include <iomanip>
#include "gauss.h"
#include <stdlib.h>
#include <iostream>

using namespace nicole;
int main()
{
	srand(time(NULL));
	const unsigned int size = 100;
	hpmat A = hprandn(10,10,size);
	hpmat b = hprandn(10,10,size,1);
	
	/*Create normal equations*/
	A = A.t()*A;
	b = A.t()*b;
	
	/*Solve*/
	std::cout << conjgrad<hp>::solve(A,b,.000001,1000);
}