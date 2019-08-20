#pragma once
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include "matrix.h"
#include "gauss.h"


namespace nicole
{
	template<class T>
	class random_generator
	{
	protected:
		long long s;
		void update_seed(void);
		long long a,c,m;
	public:
		random_generator(void);
		random_generator(const long long);
		random_generator(const long long,const long long,const long long);
		random_generator(const long long,const long long,const long long,const long long);
		random_generator(const random_generator<T>&);
		~random_generator(void);
		random_generator<T> & operator=(const random_generator<T>&);
		long long & seed(void);
		random_generator<T> & seed(const long long);
		random_generator<T> & set_generator_params(const long long,const long long,const long long);
		random_generator<T> & reset_default_params(void);
		
		/*Generations of Random Variables -> More to come :)*/
		
		/*Uniform random variables*/
		T uniform(void);
		T uniform(T,T);
		matrix<T> uniform(unsigned int);
		matrix<T> uniform(T,T,unsigned int);
		matrix<T> uniform(unsigned int,unsigned int);
		matrix<T> uniform(T,T,unsigned int,unsigned int);
		
		/*Weibull random variables*/
		T weibull(T,T);
		matrix<T> weibull(T,T,unsigned int);
		matrix<T> weibull(T,T,unsigned int,unsigned int);
		
		/*Normal random variables using Box-Muler transform*/
		T randn(void);
		matrix<T> randn(unsigned int);
		matrix<T> randn(unsigned int,unsigned int);
		T randn(T,T);
		matrix<T> randn(T,T,unsigned int);
		matrix<T> randn(T,T,unsigned int,unsigned int);
		
		/*Normal Random variables using Marsaglias Polar Method*/
		T randn_pole(void);
		matrix<T> randn_pole(unsigned int);
		matrix<T> randn_pole(unsigned int,unsigned int);
		T randn_pole(T,T);
		matrix<T> randn_pole(T,T,unsigned int);
		matrix<T> randn_pole(T,T,unsigned int,unsigned int);
		
		/*Generate from a multivaiate gaussian distrbution*/
		matrix<T> multi_gaussian(const matrix<T>&,const matrix<T>&); //Mean vector and cov matrix
		matrix<T> multi_gaussian(const matrix<T>&,const matrix<T>&,unsigned int);
		
		/*Generate Expodential RV*/
		T exponential(T);
		matrix<T> exponential(T,unsigned int);
		matrix<T> exponential(T,unsigned int,unsigned int);
		/*For stratification simulations*/
		T exponential(T,T,T);
		
		/*Generate a standard die*/
		T die(void);
		T die(double);
		matrix<T> die(unsigned int);
		matrix<T> die(unsigned int,unsigned int);
		
		/*Generate coin flip*/
		T coin(void);
		matrix<T> coin(unsigned int);
		matrix<T> coin(unsigned int,unsigned int);
	
	};
	
	template<class T>
	random_generator<T>::random_generator(void)
	{
		s = time(NULL);
		a = (long long) 1664525;
		c = (long long) 1013904223;
		m = (long long) std::pow((long long)2,(long long)32);
		if(s >= m) {std::cout << "Runtime Error\nSeed must be less then: " << m << "\n" ; exit(1);} else {}
	}
	
	template<class T>
	random_generator<T>::random_generator(const long long se)
	{
		s = se;
		a = (long long) 1664525;
		c = (long long) 1013904223;
		m = (long long) std::pow((long long)2,(long long)32);
		if(s >= m) {std::cout << "Runtime Error\nSeed must be less then: " << m << "\n"; exit(1);} else {}
	}
	
	template<class T>
	random_generator<T>::random_generator(const long long ap,const long long cp,const long long mp)
	{
		s = time(NULL);
		a = ap;
		c = cp;
		m = mp;
		if(s >= m) {std::cout << "Runtime Error\nSeed must be less then: " << m << "\n" ; exit(1);} else {}
	}
	
	template<class T>
	random_generator<T>::random_generator(const long long se,const long long ap,const long long cp,const long long mp)
	{
		s = se;
		a = ap;
		c = cp;
		m = mp;
		if(s >= m) {std::cout << "Runtime Error\nSeed must be less then: " << m << "\n" ; exit(1);} else {}
	}
	
	template<class T>
	random_generator<T>::random_generator(const random_generator<T> &rhs)
	{
		s = rhs.s;
		a = rhs.a;
		c = rhs.c;
		m = rhs.m;
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
		a = rhs.a;
		c = rhs.c;
		m = rhs.m;
		return *this;
	}
	
	template<class T>
	inline
	long long & random_generator<T>::seed(void)
	{
		return s;
	}
	
	template<class T>
	inline
	random_generator<T> & random_generator<T>::seed(const long long se)
	{
		s = se;
		return *this;
	}
	
	template<class T>
	inline
	random_generator<T> & random_generator<T>::set_generator_params(const long long ap,const long long cp,const long long mp)
	{
		a = ap;
		c = cp;
		m = mp;
		if(s >= m) {std::cout << "Runtime Error\nSeed must be less then: " << m << "\n" ; exit(1);} else {}
		return *this;
	}
	
	template<class T>
	inline
	void random_generator<T>::update_seed(void)
	{
		s = (a*s + c) % m;;
	}
	
	template<class T>
	inline
	random_generator<T> & random_generator<T>::reset_default_params(void)
	{
		a = (long long) 1664525;
		c = (long long) 1013904223;
		m = (long long) std::pow((long long)2,(long long)32);
		if(s >= m) {std::cout << "Runtime Error\nSeed must be less then: " << m << "\n" ; exit(1);} else {}
		return *this;
	}
	
	template<class T>
	inline
	T random_generator<T>::uniform(void)
	{
		T x = (T) s/m;
		update_seed();
		return x;
	}
	
	template<class T>
	inline
	T random_generator<T>::uniform(T a,T b)
	{
		return (b-a)*uniform() + a;
	}
	
	template<class T>
	matrix<T> random_generator<T>::uniform(unsigned int N)
	{
		matrix<T> X(N,N);

		for(int i = 0; i < N; i++)
		{
			for(int j = 0; j < N; j++)
			{
				X(i,j) = uniform();
			}
		}
		return X;
	}
	
	template<class T>
	matrix<T> random_generator<T>::uniform(unsigned int n,unsigned int m)
	{
		matrix<T> X(n,m);
		for(int i = 0; i < n; i++)
		{
			for(int j = 0; j < m; j++)
			{
				X(i,j) = uniform();
			}
		}
		return X;
	}
	
	template<class T>
	matrix<T> random_generator<T>::uniform(T a,T b,unsigned int n,unsigned int m)
	{
		matrix<T> X(n,m);
		for(int i = 0; i < n; i++)
		{
			for(int j = 0; j < m; j++)
			{
				X(i,j) = uniform(a,b);
			}
		}
		return X;
	}
	
	template<class T>
	inline
	T random_generator<T>::weibull(T alpha,T beta)
	{
		return std::pow(-std::log(uniform())/alpha,1./beta);
	}
	
	template<class T>
	matrix<T> random_generator<T>::weibull(T alpha,T beta,unsigned int N)
	{
		matrix<T> X(N,N);
		for(int i = 0; i < N; i++)
		{
			for(int j = 0; j < N; j++)
			{
				X(i,j) = weibull(alpha,beta);
			}
		}
		return X;
	}
	
	template<class T>
	matrix<T> random_generator<T>::weibull(T alpha,T beta,unsigned int n,unsigned int m)
	{
		matrix<T> X(n,m);
		for(int i = 0; i < n; i++)
		{
			for(int j = 0; j < m; j++)
			{
				X(i,j) = weibull(alpha,beta);
			}
		}
		return X;
	}
	
	/*Random standard normal class!*/
	template<class T>
	inline
	T random_generator<T>::randn(void)
	{
		return std::sqrt(-(T)2*std::log(uniform()))*std::cos((T)2.*(T)pi*uniform());
	}
	
	template<class T>
	matrix<T> random_generator<T>::randn(unsigned int N)
	{
		/*Initalize a matrix to store the data*/
		matrix<T> X(N,N);
		/*Get a pointer to the data*/
		T *raw_dat = X.raw_data();
		/*Generate the random variables*/
		for(int i = 0; i < N*N-1; i++)
		{
			T u1 = uniform(); T u2 = uniform(); 
			T R = std::sqrt(-(T)2*std::log(u1));
			raw_dat[i] = R*std::cos((T)2.*(T)pi*u2);
			raw_dat[i+1] = R*std::sin((T)2.*(T)pi*u2);
		}
		raw_dat[N*N-1] = std::sqrt(-(T)2*std::log(uniform()))*std::sin((T)2.*(T)pi*uniform());
		return X;
	}
	
	template<class T>
	matrix<T> random_generator<T>::randn(unsigned int n,unsigned int m)
	{
		/*Initalize a matrix to store the data*/
		matrix<T> X(n,m);
		/*Get a pointer to the data*/
		T *raw_dat = X.raw_data();
		/*Generate the random variables*/
		for(int i = 0; i < n*m-1; i++)
		{
			T u1 = uniform(); T u2 = uniform();  
			T R = std::sqrt(-(T)2*std::log(u1));
			raw_dat[i] = R*std::cos((T)2.*(T)pi*u2);
			raw_dat[i+1] = R*std::sin((T)2.*(T)pi*u2);
		}
		raw_dat[m*n-1] = std::sqrt(-(T)2*std::log(uniform()))*std::sin((T)2.*(T)pi*uniform());
		return X;
	}
	
	template<class T>
	inline
	T random_generator<T>::randn(T mu,T var)
	{
		return std::sqrt(var)*randn() + mu;
	}
	
	template<class T>
	matrix<T> random_generator<T>::randn(T mu,T var,unsigned int N)
	{
		return std::sqrt(var)*randn(N) + mu;
	}
	
	template<class T>
	matrix<T> random_generator<T>::randn(T mu,T var,unsigned int n,unsigned int m)
	{
		return std::sqrt(var)*randn(n,m) + mu;
	}
	
	template<class T>
	T random_generator<T>::randn_pole(void)
	{
		T S = (T) 2;
		T u1,u2,R;
		/*Accepts the random cordinates*/
		while(S > (T)1)
		{
			u1 = uniform(-(T)1,(T)1); u2 = uniform(-(T)1,(T)1); 
			S = (u1*u1 + u2*u2);
		}
		return u1*std::sqrt(-(T)2*std::log(S)/S);
	}
	
	template<class T>
	matrix<T> random_generator<T>::randn_pole(unsigned int N)
	{
		/*Initalize a matrix to store the data*/
		matrix<T> X(N,N);
		/*Get a pointer to the data*/
		T *raw_dat = X.raw_data();
		/*Generate the random variables*/
		for(int i = 0; i < N*N-1; i++)
		{
			T S = (T) 2;
			T u1,u2,R;
			/*Accepts the random cordinates*/
			while(S > (T)1)
			{
				u1 = uniform(-(T)1,(T)1); u2 = uniform(-(T)1,(T)1); 
				S = (u1*u1 + u2*u2);
			}
			R = std::sqrt(-(T)2*std::log(S)/S);
			raw_dat[i] = u1*R;
			raw_dat[i+1] = u2*R;
		}
		raw_dat[N*N-1] = randn_pole();
		return X;
	}
	
	template<class T>
	matrix<T> random_generator<T>::randn_pole(unsigned int n,unsigned int m)
	{
		/*Initalize a matrix to store the data*/
		matrix<T> X(n,m);
		/*Get a pointer to the data*/
		T *raw_dat = X.raw_data();
		for(int i = 0; i < n*m-1; i++)
		{
			T S = (T) 2;
			T u1,u2,R;
			/*Accepts the random cordinates*/
			while(S > (T)1)
			{
				u1 = uniform(-(T)1,(T)1); u2 = uniform(-(T)1,(T)1); 
				S = (u1*u1 + u2*u2);
			}
			R = std::sqrt(-(T)2*std::log(S)/S);
			raw_dat[i] = u1*R;
			raw_dat[i+1] = u2*R;
		}
		raw_dat[n*m-1] = randn_pole();
		return X;
	}
	
	template<class T>
	T random_generator<T>::randn_pole(T mu,T var)
	{
		return std::sqrt(var)*randn_pole() + mu;
	}
	
	template<class T>
	matrix<T> random_generator<T>::randn_pole(T mu,T var,unsigned int N)
	{
		return std::sqrt(var)*randn_pole(N) + mu;
	}
	
	template<class T>
	matrix<T> random_generator<T>::randn_pole(T mu,T var,unsigned int n,unsigned int m)
	{
		return std::sqrt(var)*randn_pole(n,m) + mu;
	}
	
	template<class T>
	matrix<T> random_generator<T>::multi_gaussian(const matrix<T> &mu,const matrix<T> &C)
	{
		if(!C.is_square() || !mu.is_col_vector() || mu.numel() != C.rows()) {std::cout << "Runtime Error\nInvalid Sizes\n"; exit(1);} else {}
		/*Initalize the X vector*/
		matrix<T> X = matrix<T>::zeros(mu.numel(),(ui)1);
		/*Generate independent gaussians*/
		matrix<T> Z = randn(mu.numel(),(ui)1);
		/*Reduce C into A*A' = C*/
		matrix<T> A = gauss<T>::cholesky(C);
		return mu + A*Z;
		
	}
	
	template<class T>
	matrix<T> random_generator<T>::multi_gaussian(const matrix<T> &mu,const matrix<T> &C,unsigned int N)
	{
		if(!C.is_square() || !mu.is_col_vector() || mu.numel() != C.rows()) {std::cout << "Runtime Error\nInvalid Sizes\n"; exit(1);} else {}
		/*Initalize the X vector*/
		matrix<T> X = matrix<T>::zeros(mu.numel(),N);
		/*Reduce C into A*A' = C*/
		matrix<T> A = gauss<T>::cholesky(C);
		
		/*Iterate and initalize all the columns*/
		for(int i = 0; i < N; i++)
		{
			/*Generate independent gaussians*/
			matrix<T> Z = randn(mu.numel(),(ui)1);
			/*Our Z variable*/
			matrix<T> V = mu + A*Z;
			/*Add to matrix*/
			for(int j = 0; j < mu.numel(); j++)
			{
				X(j,i) = V(j);
			}
		}
		return X;
	}
	
	/*Exp distribution generation*/
	template<class T>
	inline
	T random_generator<T>::exponential(T lambda)
	{
		return -std::log(uniform())/lambda;
	}
	
	/*For stratification*/
	template<class T>
	inline
	T random_generator<T>::exponential(T lambda,T a,T b)
	{
		return -std::log(uniform(a,b))/lambda;
	}
	
	template<class T>
	matrix<T> random_generator<T>::exponential(T lambda,unsigned int N)
	{
		return -log(uniform(N))/lambda;
	}
	
	template<class T>
	matrix<T> random_generator<T>::exponential(T lambda,unsigned int n,unsigned int m)
	{
		return -log(uniform(n,m))/lambda;
	}
	
	template<class T>
	T random_generator<T>::die(void)
	{
		return std::floor(uniform((T)0,(T)6)) + 1;
	}
	
	template<class T>
	T random_generator<T>::die(double u)
	{
		return std::floor((T)6*u) + (T)1;
	}
	
	template<class T>
	matrix<T> random_generator<T>::die(unsigned int N)
	{
		return floor((T)6*uniform(N)) + (T)1;
	}
	
	template<class T>
	matrix<T> random_generator<T>::die(unsigned int n,unsigned int m)
	{
		return floor((T)6*uniform(n,m)) + (T)1;
	}
	
	template<class T>
	T random_generator<T>::coin(void)
	{
		if(uniform() > (T).5) {return (T)1.;} else {return (T)0.;}
	}
	
	template<class T>
	matrix<T> random_generator<T>::coin(unsigned int N)
	{
		matrix<T> C = matrix<T>::zeros(N,N);
		for(int i = 0; i < N*N;++i)
		{
			C.raw_data()[i] = coin();
		}
		return C;
	}
	
	template<class T>
	matrix<T> random_generator<T>::coin(unsigned int n,unsigned int m)
	{
		matrix<T> C = matrix<T>::zeros(n,m);
		for(int i = 0; i < n*m;++i)
		{
			C.raw_data()[i] = coin();
		}
		return C;
	}
	
}