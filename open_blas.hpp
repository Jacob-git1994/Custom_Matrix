#pragma once
#include <cblas.h>
#include "matrix.h"
#include <iostream>
#include <stdlib.h>

/*Compile and link the the open blas folder -> dont forget to use make on it - Documentation somewhere online*/
/*-I/Users/Jacob/Documents/Math/openBLAS/OpenBLAS-0.2.20 -L//Users/Jacob/Documents/Math/openBLAS/OpenBLAS-0.2.20 -lopenblas -pthread -lgfortran*/

namespace nicole
{
	template<class T>
	class cblas
	{
	public:
		static matrix<T> cblas_gemm(const matrix<T>&,const matrix<T>&,const matrix<T>&,T,T);
	};
	
	template<>
	matrix<double> cblas<double>::cblas_gemm(const matrix<double> &A,const matrix<double> &B,const matrix<double> &C,double alpha,double beta)
	{
		/*Cheak to see if the method will work*/
		if(A.cols() != B.rows() || A.rows() != C.rows() || B.cols() != C.cols()) {std::cout << "Runtime Error\nMatrix dimention mismatch\n"; exit(1);}
		else {}
		
		/*Initalize parameters*/
		const unsigned int M = A.rows();
		const unsigned int N = B.cols();
		const unsigned int K = A.cols();
		
		/*Copy the Matrix*/
		const matrix<double> A_C = A;
		const matrix<double> B_C = B;
		const matrix<double> C_C = C;
		
		/*Multiply the system*/
		cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,M,N,K,alpha,A_C.raw_data(),K,B_C.raw_data(),N,beta,C_C.raw_data(),N);
		return C_C;
	}

	template<>
	matrix<float> cblas<float>::cblas_gemm(const matrix<float> &A,const matrix<float> &B,const matrix<float> &C,float alpha,float beta)
	{
		/*Cheak to see if the method will work*/
		if(A.cols() != B.rows() || A.rows() != C.rows() || B.cols() != C.cols()) {std::cout << "Runtime Error\nMatrix dimention mismatch\n"; exit(1);}
		else {}
		
		/*Initalize parameters*/
		const unsigned int M = A.rows();
		const unsigned int N = B.cols();
		const unsigned int K = A.cols();
		
		/*Copy the Matrix*/
		const matrix<float> A_C = A;
		const matrix<float> B_C = B;
		const matrix<float> C_C = C;
		
		/*Multiply the system*/
		cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,M,N,K,alpha,A_C.raw_data(),K,B_C.raw_data(),N,beta,C_C.raw_data(),N);
		return C_C;
	}












}