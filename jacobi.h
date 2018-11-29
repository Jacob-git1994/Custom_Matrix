#pragma once
#include "matrix.h"
#include <iostream>
#include <stdlib.h>

namespace nicole
{
	/*Simple implementation of the jacobi method for use with my matrix class*/
	template<class T>
	class jacobi
	{
	public:
		static matrix<T> solve(const matrix<T>&,const matrix<T>&);
		static matrix<T> solve(const matrix<T>&,const matrix<T>&,T);
		static matrix<T> solve(const matrix<T>&,const matrix<T>&,unsigned int);
		static matrix<T> solve(const matrix<T>&,const matrix<T>&,T,unsigned int);
	};
	
	template<class T>
	matrix<T> jacobi<T>::solve(const matrix<T> &A,const matrix<T> &b)
	{
		if(A.rows() != b.rows()) {std::cout << "Runtime Error\nMatrix not the same size\n"; exit(1);} else {}
		/*Initalize default jacobi paremeters*/
		const unsigned int max_iter = 100;
		const T tol = .001;
		/*Initalize our guess*/
		matrix<T> g1 = b;
		matrix<T> g2 = matrix<T>::zeros(b.numel(),1);
		/*Initalize a counter*/
		unsigned int step = 0;
		/*Start solving*/
		while(norm_2(g1 - g2)/norm_2(g1) > tol && step < max_iter)
		{
			g1 = g2;
#pragma omp parallel for simd
			for(int i = 0; i < b.numel(); i++)
			{
				T sigma = 0.;
				for(int j = 0; j < b.numel(); j++)
				{
					if(i != j) {sigma += A(i,j)*g2(j);}
					else if(i == j) {continue;}
					else {std::cout << "Runtime Error\nSomething went wrong\n";exit(1);}
				}
				g2(i) = (b(i) - sigma)/A(i,i);
			}
			step++;
		}
		std::cout << "Solved with " << step << " iterations\n";
		return g2;
	}
	
	template<class T>
	matrix<T> jacobi<T>::solve(const matrix<T> &A,const matrix<T> &b,T tol)
	{
		if(A.rows() != b.rows()) {std::cout << "Runtime Error\nMatrix not the same size\n"; exit(1);} else {}
		/*Initalize default jacobi paremeters*/
		const unsigned int max_iter = 100;
		/*Initalize our guess*/
		matrix<T> g1 = b;
		matrix<T> g2 = matrix<T>::zeros(b.numel(),1);
		/*Initalize a counter*/
		unsigned int step = 0;
		/*Start solving*/
		while(norm_2(g1 - g2)/norm_2(g1) > tol && step < max_iter)
		{
			g1 = g2;
#pragma omp parallel for simd
			for(int i = 0; i < b.numel(); i++)
			{
				T sigma = 0.;
				for(int j = 0; j < b.numel(); j++)
				{
					if(i != j) {sigma += A(i,j)*g2(j);}
					else if(i == j) {continue;}
					else {std::cout << "Runtime Error\nSomething went wrong\n";exit(1);}
				}
				g2(i) = (b(i) - sigma)/A(i,i);
			}
			step++;
		}
		std::cout << "Solved with " << step << " iterations\n";
		return g2;
	}
	
	template<class T>
	matrix<T> jacobi<T>::solve(const matrix<T> &A,const matrix<T> &b,unsigned int max_iter)
	{
		if(A.rows() != b.rows()) {std::cout << "Runtime Error\nMatrix not the same size\n"; exit(1);} else {}
		/*Initalize default jacobi paremeters*/
		const T tol = .001;
		/*Initalize our guess*/
		matrix<T> g1 = b;
		matrix<T> g2 = matrix<T>::zeros(b.numel(),1);
		/*Initalize a counter*/
		unsigned int step = 0;
		/*Start solving*/
		while(norm_2(g1 - g2)/norm_2(g1) > tol && step < max_iter)
		{
			g1 = g2;
#pragma omp parallel for simd
			for(int i = 0; i < b.numel(); i++)
			{
				T sigma = 0.;
				for(int j = 0; j < b.numel(); j++)
				{
					if(i != j) {sigma += A(i,j)*g2(j);}
					else if(i == j) {continue;}
					else {std::cout << "Runtime Error\nSomething went wrong\n";exit(1);}
				}
				g2(i) = (b(i) - sigma)/A(i,i);
			}
			step++;
		}
		std::cout << "Solved with " << step << " iterations\n";
		return g2;
	}
	
	template<class T>
	matrix<T> jacobi<T>::solve(const matrix<T> &A,const matrix<T> &b,T tol,unsigned int max_iter)
	{
		if(A.rows() != b.rows()) {std::cout << "Runtime Error\nMatrix not the same size\n"; exit(1);} else {}
		/*Initalize our guess*/
		matrix<T> g1 = b;
		matrix<T> g2 = matrix<T>::zeros(b.numel(),1);
		/*Initalize a counter*/
		unsigned int step = 0;
		/*Start solving*/
		while(norm_2(g1 - g2)/norm_2(g1) > tol && step < max_iter)
		{
			g1 = g2;
#pragma omp parallel for simd
			for(int i = 0; i < b.numel(); i++)
			{
				T sigma = 0.;
				for(int j = 0; j < b.numel(); j++)
				{
					if(i != j) {sigma += A(i,j)*g2(j);}
					else if(i == j) {continue;}
					else {std::cout << "Runtime Error\nSomething went wrong\n";exit(1);}
				}
				g2(i) = (b(i) - sigma)/A(i,i);
			}
			step++;
		}
		std::cout << "Solved with " << step << " iterations\n";
		return g2;
	}
}