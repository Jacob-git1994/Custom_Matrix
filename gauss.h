#pragma once
#include <iostream>
#include <stdlib.h>
#include <cmath>
#include "matrix.h"

namespace nicole
{
	//Checks to see if we want a permuation
	enum piviot_args {nopiviot = 0,piviot = 1};

	template<class T>
	class gauss
	{
	public:
		/*Gaussian solver args -- make sure the piviot vector is uninitalized*/
		static matrix<T> solve(const matrix<T>&,const matrix<T>&);
		static matrix<T> solve(const matrix<T>&,const matrix<T>&,piviot_args);
		static matrix<T> solve(const matrix<T>&,const matrix<T>&,matrix<unsigned int>&);
		static matrix<T> solve(const matrix<T>&,const matrix<T>&,matrix<unsigned int>&,piviot_args);
		
		/*Banded Matrix Solver*/
		static matrix<T> solve_banded(const matrix<T>&,const matrix<T>&,const unsigned int);

		/*Using gaussian elemination returns the determinate*/
		static T det(const matrix<T>&);
	};

	/*Basic solver with piviot on*/
	template<class T>
	matrix<T> gauss<T>::solve(const matrix<T> &A,const matrix<T> &b)
	{
		/*Check to make sure sizes*/
		if(A.cols() != b.rows() || !A.is_square()) {std::cout << "Runtime Error\nMatricies not the same size\n"; exit(1);} else {}
		
		/*Initalize copies and generate a vector for the solution*/
		matrix<T> A_temp = A;
		matrix<T> b_temp = b;
		matrix<T> sol(b.numel(),1);
		
		/*Start solving*/
		for(int k = 0; k < A_temp.cols(); k++)
		{
			/*Generate where to pivot from*/
			unsigned int max_index = 0;
			T max = (T) 0.0;
#pragma omp parallel for simd
			for(int i = k; i < A_temp.rows(); i++)
			{
				if(std::fabs(A_temp(i,k)) > max) {max_index = i; max = std::fabs(A_temp(i,k));}
				else {}
			}
			if(max_index != k)
			{
				A_temp.row_swap(k,max_index);
				b_temp.row_swap(k,max_index);
			}
			
			/*Start parallel solving*/
#pragma omp parallel for simd
			for(int i = k+1; i < A_temp.rows(); i++)
			{
				T aik = A_temp(i,k);
				b_temp(i) -= (b_temp(k)*aik)/A_temp(k,k);
				for(int j = 0; j < A_temp.cols(); j++)
				{
					A_temp(i,j) -= (A_temp(k,j))*aik/A_temp(k,k);
				}
			}
		}
		/*now time for backsubsitution*/
		for(int i = A_temp.rows()-1; i >= 0; i--)
		{
#pragma omp for simd ordered
			for(int j = i+1; j < A_temp.rows(); j++)
			{
				b_temp(i) -= sol(j)*A_temp(i,j); 
			}
			sol(i) = b_temp(i)/A_temp(i,i);
		}
		return sol;
	}

	template<class T>
	matrix<T> gauss<T>::solve(const matrix<T> &A,const matrix<T> &b,piviot_args pv)
	{
		/*Check to make sure sizes*/
		if(A.cols() != b.rows() || !A.is_square()) {std::cout << "Runtime Error\nMatricies not the same size\n"; exit(1);} else {}
		
		/*Initalize copies and generate a vector for the solution*/
		matrix<T> A_temp = A;
		matrix<T> b_temp = b;
		matrix<T> sol(b.numel(),1);
		
		/*Start solving*/
		for(int k = 0; k < A_temp.cols(); k++)
		{
			if(pv)
			{
				/*Generate where to pivot from*/
				unsigned int max_index = 0;
				T max = (T) 0.0;
#pragma omp parallel for simd
				for(int i = k; i < A_temp.rows(); i++)
				{
					if(std::fabs(A_temp(i,k)) > max) {max_index = i; max = std::fabs(A_temp(i,k));}
					else {}
				}
				if(max_index != k)
				{
					A_temp.row_swap(k,max_index);
					b_temp.row_swap(k,max_index);
				}
			else if(pv)
			{
					/*Do nothing*/
			}
			else
			{
				std::cout << "Runtime Error\nInvlaid piviot args\n"; exit(1);
			}
			}
			
			/*Start parallel solving*/
#pragma omp parallel for simd
			for(int i = k+1; i < A_temp.rows(); i++)
			{
				T aik = A_temp(i,k);
				b_temp(i) -= (b_temp(k)*aik)/A_temp(k,k);
				for(int j = 0; j < A_temp.cols(); j++)
				{
					A_temp(i,j) -= (A_temp(k,j))*aik/A_temp(k,k);
				}
			}
		}
		/*now time for backsubsitution*/
		for(int i = A_temp.rows()-1; i >= 0; i--)
		{
#pragma omp for simd ordered
			for(int j = i+1; j < A_temp.rows(); j++)
			{
				b_temp(i) -= sol(j)*A_temp(i,j); 
			}
			sol(i) = b_temp(i)/A_temp(i,i);
		}
		return sol;
	}

	template<class T>
	matrix<T> gauss<T>::solve(const matrix<T> &A,const matrix<T> &b,matrix<unsigned int> &perm)
	{
		/*Check to make sure sizes*/
		if(A.cols() != b.rows() || !A.is_square()) {std::cout << "Runtime Error\nMatricies not the same size\n"; exit(1);} else {}
		
		/*Initalize copies and generate a vector for the solution*/
		matrix<T> A_temp = A;
		matrix<T> b_temp = b;
		matrix<T> sol(b.numel(),1);
		
		/*Initalizes the perm vector*/
		perm = nicole::matrix<unsigned int>::linspace(0,A.rows()-1,A.rows());

		/*Start solving*/
		for(int k = 0; k < A_temp.cols(); k++)
		{
			/*Generate where to pivot from*/
			unsigned int max_index = 0;
			T max = (T) 0.0;
#pragma omp parallel for simd
			for(int i = k; i < A_temp.rows(); i++)
			{
				if(std::fabs(A_temp(i,k)) > max) {max_index = i; max = std::fabs(A_temp(i,k));}
				else {}
			}
			if(max_index != k)
			{
				A_temp.row_swap(k,max_index);
				b_temp.row_swap(k,max_index);
				perm.col_swap(k,max_index);
			}
			
			/*Start parallel solving*/
#pragma omp parallel for simd
			for(int i = k+1; i < A_temp.rows(); i++)
			{
				T aik = A_temp(i,k);
				b_temp(i) -= (b_temp(k)*aik)/A_temp(k,k);
				for(int j = 0; j < A_temp.cols(); j++)
				{
					A_temp(i,j) -= (A_temp(k,j))*aik/A_temp(k,k);
				}
			}
		}
		/*now time for backsubsitution*/
		for(int i = A_temp.rows()-1; i >= 0; i--)
		{
#pragma omp for simd ordered
			for(int j = i+1; j < A_temp.rows(); j++)
			{
				b_temp(i) -= sol(j)*A_temp(i,j); 
			}
			sol(i) = b_temp(i)/A_temp(i,i);
		}
		return sol;
	}

	template<class T>
	matrix<T> gauss<T>::solve(const matrix<T> &A,const matrix<T> &b,matrix<unsigned int> &perm,piviot_args pv)
	{
		/*Check to make sure sizes*/
		if(A.cols() != b.rows() || !A.is_square()) {std::cout << "Runtime Error\nMatricies not the same size\n"; exit(1);} else {}
		
		/*Initalize copies and generate a vector for the solution*/
		matrix<T> A_temp = A;
		matrix<T> b_temp = b;
		matrix<T> sol(b.numel(),1);
		
		/*Initalizes the perm vector*/
		perm = nicole::matrix<unsigned int>::linspace(0,A.rows()-1,A.rows());
		
		/*Start solving*/
		for(int k = 0; k < A_temp.cols(); k++)
		{
			if(pv)
			{
				/*Generate where to pivot from*/
				unsigned int max_index = 0;
				T max = (T) 0.0;
#pragma omp parallel for simd
				for(int i = k; i < A_temp.rows(); i++)
				{
					if(std::fabs(A_temp(i,k)) > max) {max_index = i; max = std::fabs(A_temp(i,k));}
					else {}
				}
				if(max_index != k)
				{
					A_temp.row_swap(k,max_index);
					b_temp.row_swap(k,max_index);
					perm.col_swap(k,max_index);
				}
			else if(pv)
			{
					/*Do nothing*/
			}
			else
			{
				std::cout << "Runtime Error\nInvlaid piviot args\n"; exit(1);
			}
			}
			
			/*Start parallel solving*/
#pragma omp parallel for simd
			for(int i = k+1; i < A_temp.rows(); i++)
			{
				T aik = A_temp(i,k);
				b_temp(i) -= (b_temp(k)*aik)/A_temp(k,k);
				for(int j = 0; j < A_temp.cols(); j++)
				{
					A_temp(i,j) -= (A_temp(k,j))*aik/A_temp(k,k);
				}
			}
		}
		/*now time for backsubsitution*/
		for(int i = A_temp.rows()-1; i >= 0; i--)
		{
#pragma omp for simd ordered
			for(int j = i+1; j < A_temp.rows(); j++)
			{
				b_temp(i) -= sol(j)*A_temp(i,j); 
			}
			sol(i) = b_temp(i)/A_temp(i,i);
		}
		return sol;
	}
		
	template<class T>
	T gauss<T>::det(const matrix<T> &A)
	{
		/*Check to make sure sizes*/
		if(!A.is_square()) {std::cout << "Runtime Error\nMatrix must be square\n"; exit(1);} else {}
		
		/*Determinate tester*/
		double alpha = 1;
		
		/*Initalize copy of A*/
		matrix<T> A_temp = A;
		
		/*Start solving*/
		for(int k = 0; k < A_temp.cols(); k++)
		{
			/*Generate where to pivot from*/
			unsigned int max_index = 0;
			T max = (T) 0.0;
#pragma omp parallel for simd
			for(int i = k; i < A_temp.rows(); i++)
			{
				if(std::fabs(A_temp(i,k)) > max) {max_index = i; max = std::fabs(A_temp(i,k));}
				else {}
			}
			if(max_index != k)
			{
				A_temp.row_swap(k,max_index);
				alpha*=-1;
			}
			
			/*Start parallel solving*/
#pragma omp parallel for simd
			for(int i = k+1; i < A_temp.rows(); i++)
			{
				T aik = A_temp(i,k);
				for(int j = 0; j < A_temp.cols(); j++)
				{
					A_temp(i,j) -= (A_temp(k,j))*aik/A_temp(k,k);
				}
			}
		}
		/*Multiply down the diagonal*/
		T det = (T) 1.0;
#pragma omp parallel for simd reduction(*:det)
		for(int i = 0; i < A_temp.rows(); i++)
		{
			det *= A_temp(i,i);
		}
		return alpha*det;
	}
}