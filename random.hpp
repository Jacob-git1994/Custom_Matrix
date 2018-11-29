#pragma once
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include "matrix.h"


namespace nicole
{
	template<class T>
	class random_generator
	{
	private:
		long long s;
		void update_seed(void);
	public:
		random_generator(void);
		random_generator(const long long);
		~random_generator(void);
		random_generator<T> & operator=(const random_generator<T>&);
		long long & seed(void);
		
		/*Generations of Random Variables*/
		/*Uniform random variables*/
		T uniform(void);
		T uniform(T,T);
		matrix<T> uniform(unsigned int);
		matrix<T> uniform(T,T,unsigned int);
		matrix<T> uniform(unsigned int,unsigned int);
		matrix<T> uniform(T,T,unsigned int,unsigned int);
	
	};
	
	template<class T>
	random_generator<T>::random_generator(void)
	{
		s = time(NULL);
	}
	
	template<class T>
	random_generator<T>::random_generator(const long long se)
	{
		if(se >= (long long)std::pow(2,32)) {std::cout << "Runtime Error\nSeed must be less then 2^32\n"; exit(1);} else {}
		s = se;
	}
	
	template<class T>
	random_generator<T>::~random_generator(void)
	{
		/*Do nothing*/
	}
	
	template<class T>
	random_generator<T> & random_generator<T>::operator=(const random_generator<T> &rhs)
	{
		s = rhs.s;
		return *this;
	}
	
	template<class T>
	long long & random_generator<T>::seed(void)
	{
		return s;
	}
	
	template<class T>
	void random_generator<T>::update_seed(void)
	{
		s = ((long long)1664525*s + (long long)1013904223) % (long long)std::pow(2,32);
	}
	
	template<class T>
	T random_generator<T>::uniform(void)
	{
		T x =  (T)s/(T)std::pow(2,32);
		update_seed();
		return x;
	}
	
	template<class T>
	T random_generator<T>::uniform(T a,T b)
	{
		return (b-a)*uniform() + a;
	}
	
	template<class T>
	matrix<T> random_generator<T>::uniform(unsigned int n,unsigned int m)
	{
		matrix<T> X = matrix<T>::zeros(n,m);
#pragma omp parallel for simd
		for(int i = 0; i < n; i++)
		{
			for(int j = 0; j < m; j++)
			{
				X(i,j) = uniform();
			}
		}
		return X;
	}
	
	
	
	
	
}