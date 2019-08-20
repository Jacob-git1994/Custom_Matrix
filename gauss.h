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
		static matrix<T> solve_tridiag(const matrix<T>&,const matrix<T>&);
		
		/*Choleski-Factorization is a modifued gaussian elimination*/
		static matrix<T> cholesky(const matrix<T>&); // Returns L where L*L' = A

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
#pragma omp parallel for
			for(int i = k+1; i < A_temp.rows(); i++)
			{
				T aik = A_temp(i,k);
				b_temp(i) -= (b_temp(k)*aik)/A_temp(k,k);
				for(int j = k+1; j < A_temp.cols(); j++)
				{
					A_temp(i,j) -= (A_temp(k,j))*aik/A_temp(k,k);
				}
			}
		}
		/*now time for backsubsitution*/
		for(int i = A_temp.rows()-1; i >= 0; i--)
		{
#pragma omp for ordered
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
#pragma omp parallel for
			for(int i = k+1; i < A_temp.rows(); i++)
			{
				T aik = A_temp(i,k);
				b_temp(i) -= (b_temp(k)*aik)/A_temp(k,k);
				for(int j = k+1; j < A_temp.cols(); j++)
				{
					A_temp(i,j) -= (A_temp(k,j))*aik/A_temp(k,k);
				}
			}
		}
		/*now time for backsubsitution*/
		for(int i = A_temp.rows()-1; i >= 0; i--)
		{
#pragma omp for ordered
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
#pragma omp parallel for 
			for(int i = k+1; i < A_temp.rows(); i++)
			{
				T aik = A_temp(i,k);
				b_temp(i) -= (b_temp(k)*aik)/A_temp(k,k);
				for(int j = k+1; j < A_temp.cols(); j++)
				{
					A_temp(i,j) -= (A_temp(k,j))*aik/A_temp(k,k);
				}
			}
		}
		/*now time for backsubsitution*/
		for(int i = A_temp.rows()-1; i >= 0; i--)
		{
#pragma omp for ordered
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
#pragma omp parallel for 
			for(int i = k+1; i < A_temp.rows(); i++)
			{
				T aik = A_temp(i,k);
				b_temp(i) -= (b_temp(k)*aik)/A_temp(k,k);
				for(int j = k+1; j < A_temp.cols(); j++)
				{
					A_temp(i,j) -= (A_temp(k,j))*aik/A_temp(k,k);
				}
			}
		}
		/*now time for backsubsitution*/
		for(int i = A_temp.rows()-1; i >= 0; i--)
		{
#pragma omp for ordered
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
#pragma omp parallel for
			for(int i = k+1; i < A_temp.rows(); i++)
			{
				T aik = A_temp(i,k);
				for(int j = k+1; j < A_temp.cols(); j++)
				{
					A_temp(i,j) -= (A_temp(k,j))*aik/A_temp(k,k);
				}
			}
		}
		/*Multiply down the diagonal*/
		T det = (T) 1.0;
#pragma omp parallel for reduction(*:det)
		for(int i = 0; i < A_temp.rows(); i++)
		{
			det *= A_temp(i,i);
		}
		return alpha*det;
	}
	
	template<class T>
	matrix<T> gauss<T>::solve_tridiag(const matrix<T> &A,const matrix<T> &r)
	{
		/*Check to make sure sizes*/
		if(A.cols() != r.rows() || !A.is_square()) {std::cout << "Runtime Error\nMatricies not the same size\n"; exit(1);} else {}
		
		/*Break up the matrix into off diagonals and diagonals*/
		matrix<T> a = matrix<T>::zeros(A.rows(),1);
		matrix<T> b = matrix<T>::zeros(A.rows(),1);
		matrix<T> c = matrix<T>::zeros(A.rows(),1);
		matrix<T> b_c = r;
		matrix<T> sol = matrix<T>::zeros(A.rows(),1);
		
		/*Copy diag elements into vectors*/
		b(0) = A(0,0); b(A.rows()-1) = A(A.rows()-1,A.rows()-1);
		a(A.rows()-1) = A(A.rows()-1,A.rows()-2);
		c(A.rows()-1) = A(A.rows()-2,A.rows()-1); c(0) = A(0,1);
		for(int i = 1; i < A.rows()-1; i++)
		{
			a(i) = A(i,i-1);
			b(i) = A(i,i);
			c(i) = A(i,i+1);
		}
		
		/*Initalize the begining*/
		c(0) = c(0)/b(0);
		b_c(0) = b_c(0)/b(0);
		
		/*Iterate through the arrays*/
#pragma omp parallel for
		for(int i = 1; i < A.rows(); i++)
		{
			T id = (T)1.0/(b(i) - c(i-1)*a(i));
			c(i) = c(i)*id;
			b_c(i) = (b_c(i) - a(i)*b_c(i-1))*id;
		}
		
		/*Initalize the back subsitution*/
		sol(A.rows()-1) = b_c(A.rows()-1);
		/*Backsubsitution*/
		for(int i = A.rows()-2; i != -1; i--)
		{
			sol(i) = b_c(i) - c(i)*sol(i+1);
		}
		
		return sol;
	}
	
	/*Only works if matrix is symmetric pos def*/
	template<class T>
	matrix<T> gauss<T>::cholesky(const matrix<T> &A)
	{
		if(!A.is_square()) {std::cout << "Runtime Error\nA L Lp must be square\n"; exit(1);} else {}
		/*Copy A*/
		matrix<T> B = A; 
		/*Start solving*/
#pragma omp parallel for 
		for(int k = 0; k < A.rows(); k++)
		{
			B(k,k) = std::sqrt(B(k,k));
			for(int i = k+1; i < A.rows(); i++)
			{
				B(i,k) = B(i,k)/B(k,k);
			}
			for(int i = k+1; i < A.rows(); i++)
			{
				for(int j = i; j < A.rows(); j++)
				{
					B(j,i) -= B(i,k)*B(j,k);
				}
			}
		}
		/*Clean up by removing unessisary values*/
#pragma omp parallel for
		for(int i = 0; i < A.rows(); i++)
		{
			for(int j = i+1; j < A.cols(); j++)
			{
				B(i,j) = (T)0;
			}
		}
		return B;
	}
}