#pragma once
#include "matrix.h"
#include "gauss.h"
#include <iostream>
/*This class performs forming the jacobian and n-dim root finding*/

namespace nicole
{

	template<class T>
	class optimization
	{
	private:
		T step_size,tol;
		unsigned int max_iter;
		matrix<T> jacobian(matrix<T> (*)(const matrix<T>&),const matrix<T>&);
		matrix<T> inital_guess,solution,J;
	public:
		optimization(void);
		optimization(T,T);
		optimization(T,T,unsigned int);
		optimization(T,T,unsigned int,const matrix<T>&);
		~optimization(void);
		optimization<T> & set_tol(T);
		optimization<T> & set_step_size(T);
		optimization<T> & set_max_iter(unsigned int);
		optimization<T> & set_inital_guess(const matrix<T>&);
		matrix<T> get_jacobian(void);
		matrix<T> solve(matrix<T> (*)(const matrix<T>&));
		matrix<T> solve_print_stats(matrix<T> (*)(const matrix<T>&));
		
		/*Methods that allow time or added functions*/
		matrix<T> jacobian(matrix<T> (*)(const matrix<T>&,T),const matrix<T>&,T);
		matrix<T> solve(matrix<T> (*)(const matrix<T>&,T),T);
		matrix<T> solve_print_stats(matrix<T> (*)(const matrix<T>&,T),T);
	};

	template<class T>
	optimization<T>::optimization(void)
	{
		step_size = (T).001;
		tol = (T) .001;
		max_iter = 10;
	}

	template<class T>
	optimization<T>::optimization(T s,T t,unsigned int i)
	{
		step_size = s;
		tol = t;
		max_iter = i;
	}
	
	template<class T>
	optimization<T>::optimization(T s,T t,unsigned int i,const matrix<T> &ig)
	{
		step_size = s;
		tol = t;
		max_iter = i;
		if(ig.cols() != 1) {std::cout << "Runtime Error\ninital guess is not a vector\n"; exit(1);} else {}
		inital_guess = ig;
	}

	template<class T>
	optimization<T>::~optimization(void)
	{
		/*Matrix class automatically deletes*/
	}
	
	template<class T>
	optimization<T> & optimization<T>::set_tol(T t)
	{
		tol = t;
		return *this;
	}
	
	template<class T>
	optimization<T> & optimization<T>::set_step_size(T s)
	{
		step_size = s;
		return *this;
	}
	
	template<class T>
	optimization<T> & optimization<T>::set_max_iter(unsigned int i)
	{
		max_iter = i;
		return *this;
	}
	
	template<class T>
	optimization<T> & optimization<T>::set_inital_guess(const matrix<T> &inital)
	{
		if(inital.cols() != 1) {std::cout << "Runtime Error\ninital guess is not a vector\n"; exit(1);} else {}
		inital_guess = inital;
		return *this;
	}
	
	/*Jacobian Formation*/
	template<class T>
	matrix<T> optimization<T>::jacobian(matrix<T> (*func)(const matrix<T>&),const matrix<T> &g1)
	{
		const int N = inital_guess.numel();
		/*Initalze the Jacobian*/
#pragma omp parallel for
		for(int j = 0; j < N; j++)
		{
			matrix<T> copy_inital = g1;
			copy_inital(j) += step_size;
			
			if(j == 0)
			{
				J = (func(copy_inital) - func(g1))/step_size;
			}
			else if(j != 0)
			{
				J = col_concat(J,(func(copy_inital) - func(g1))/step_size);
			}
			else 
			{
				std::cout << "Runtime Error\nSomething went wrong\n"; exit(1);
			}
		}
		return J;
	}
	
	template<class T>
	matrix<T> optimization<T>::get_jacobian(void)
	{
		return J;
	}
	
	/*finding the minimum/roots*/
	template<class T>
	matrix<T> optimization<T>::solve(matrix<T> (*func)(const matrix<T>&))
	{
		matrix<T> g1 = inital_guess;
		matrix<T> g2 = inital_guess + matrix<T>::ones(inital_guess.numel(),1);
		
		unsigned int count = 0;
		while(norm_2(g2 - g1)/norm_2(g2) > tol && count++ < max_iter)
		{
			g1 = g2;
			g2 = g1 - gauss<T>::solve(jacobian(func,g1),func(g1));
		}
		return g2;
	}
	
	/*finding the minimum/roots*/
	template<class T>
	matrix<T> optimization<T>::solve_print_stats(matrix<T> (*func)(const matrix<T>&))
	{
		matrix<T> g1 = inital_guess;
		matrix<T> g2 = inital_guess + matrix<T>::ones(inital_guess.numel(),1);
		
		unsigned int count = 0;
		while(norm_2(g2 - g1)/norm_2(g2) > tol && count++ < max_iter)
		{
			g1 = g2;
			g2 = g1 - gauss<T>::solve(jacobian(func,g1),func(g1));
		}
		if(norm_2(g1 - g2)/norm_2(g2) > tol) {std::cout << "Could not find solution in tolerance and current iterations\n";}
		std::cout << "Reletive Error: " << norm_2(g1 - g2)/norm_2(g2) << "\n" << "Number of steps: " << count << "\n";
		return g2;
	}
	
	/*For matricies with extra parameter*/
	/*Jacobian Formation*/
	template<class T>
	matrix<T> optimization<T>::jacobian(matrix<T> (*func)(const matrix<T>&,T),const matrix<T> &g1,T alpha)
	{
		const int N = inital_guess.numel();
		/*Initalze the Jacobian*/
#pragma omp parallel for
		for(int j = 0; j < N; j++)
		{
			matrix<T> copy_inital = g1;
			copy_inital(j) += step_size;
			
			if(j == 0)
			{
				J = (func(copy_inital,alpha) - func(g1,alpha))/step_size;
			}
			else if(j != 0)
			{
				J = col_concat(J,(func(copy_inital,alpha) - func(g1,alpha))/step_size);
			}
			else 
			{
				std::cout << "Runtime Error\nSomething went wrong\n"; exit(1);
			}
		}
		return J;
	}
	
	/*finding the minimum/roots*/
	template<class T>
	matrix<T> optimization<T>::solve(matrix<T> (*func)(const matrix<T>&,T),T alpha)
	{
		matrix<T> g1 = inital_guess;
		matrix<T> g2 = inital_guess + matrix<T>::ones(inital_guess.numel(),1);
		
		unsigned int count = 0;
		while(norm_2(g2 - g1)/norm_2(g2) > tol && count++ < max_iter)
		{
			g1 = g2;
			g2 = g1 - gauss<T>::solve(jacobian(func,g1,alpha),func(g1,alpha));
		}
		return g2;
	}
	
	template<class T>
	matrix<T> optimization<T>::solve_print_stats(matrix<T> (*func)(const matrix<T>&,T),T alpha)
	{
		matrix<T> g1 = inital_guess;
		matrix<T> g2 = inital_guess + matrix<T>::ones(inital_guess.numel(),1);
		
		unsigned int count = 0;
		while(norm_2(g2 - g1)/norm_2(g2) > tol && count++ < max_iter)
		{
			g1 = g2;
			g2 = g1 - gauss<T>::solve(jacobian(func,g1,alpha),func(g1,alpha));
		}
		if(norm_2(g1 - g2)/norm_2(g2) > tol) {std::cout << "Could not find solution in tolerance and current iterations\n";}
		std::cout << "Reletive Error: " << norm_2(g1 - g2)/norm_2(g2) << "\n" << "Number of steps: " << count << "\n";
		return g2;
	}

}