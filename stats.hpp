#pragma once 
#include "matrix.h"
#include "random.hpp"
#include "gauss.h"
#include <stdlib.h>
#include <cmath>
#include <iomanip>

namespace nicole
{
	template<class T>
	class stats
	{
	private:
		/*Typical stats*/
		T m,v,sd;
		/*Number of times sampled for the update*/
		unsigned int cp;
		/*Random Number Generator*/
		random_generator<T> ran;
	public:
		/*Constructors*/
		stats(void);
		stats(long long);
		stats(const T,const T);
		stats(const T,const T,long long);
		stats(const stats<T>&);
		~stats(void);
		
		/*Stats Assignment operator*/
		stats<T> & operator=(const stats<T>&);
		
		/*Set the mean/Var/Count*/
		stats<T> & mean(const T);
		stats<T> & var(const T);
		stats<T> & count(const unsigned int);
		stats<T> & stdd(const T);
		
		/*Refs*/
		T & mean(void);
		T & var(void);
		T & stdd(void);
		random_generator<T> & random(void);
		unsigned int & count(void);
		
		/*Update Functions*/
		stats<T> & update(const T);
		stats<T> & reset(void);
		
		/*Static Class Stat functions ADD IN LATER*/
		/*
		static matrix<T> mean(const matrix<T>&);
		static matrix<T> var(const matrix<T>&);
		static matrix<T> stdd(const matrix<T>&);
		static matrix<T> cov(const matrix<T>&,const matrix<T>&);
		static matrix<T> cov(const matrix<T>&);
		*/
		
		/*Confidance Level*/
		T confidance_interval(const unsigned int); //Returns interval length
		/*Monte Carlo Simulations*/
		T expected(T (*)(random_generator<T>&),const long long,const unsigned int); //Fix number of realizations
		T expected(T (*)(random_generator<T>&),const T,const unsigned int); //Fix tolerance of our simulation
		
		/*Variance and mean with MSE*/
		T var(T (*)(random_generator<T>&),const long long,const long long boot,const int); //Fix number of realizations
		T expected_bootstrap(T (*)(random_generator<T>&),const long long,const long long boot,const int); //Bootstrap of mean
		
		/*Probablity estimators*/
		T greater(T (*)(random_generator<T>&),const long long,const unsigned int,const T);
		T less(T (*)(random_generator<T>&),const long long,const unsigned int,const T);
		T equal(T (*)(random_generator<T>&),const long long,const unsigned int,const T);
		T greater_equal(T (*)(random_generator<T>&),const long long,const unsigned int,const T);
		T less_equal(T (*)(random_generator<T>&),const long long,const unsigned int,const T);
		
		/*Varience reduction methods for the expected value*/
		T expected(T (*)(random_generator<T>&),const long long,const unsigned int,T (*)(random_generator<T>&),const T);
		T expected(T (*)(random_generator<T>&),const long long,const unsigned int,T (*)(random_generator<T>&),const T,T (*)(random_generator<T>&),const T);
		
		/*Calculate the expected value with strata*/
		T expected_strata(T (*)(random_generator<T>&,T,T),const long long,const unsigned int,const unsigned int);
		
		/*Common PDFS*/
		T gaussian(const T);
		T gaussian(const T,const T,const T);
		T exponential(const T,const T);
		T uniform(void);
		T uniform(const T,const T);
		
		
	};
	
	template<class T>
	stats<T>::stats(void)
	{
		m = (T)0;
		v = (T)0;
		sd = (T)0;
		cp = (ui)0;
	}
	
	template<class T>
	stats<T>::stats(const T me,const T va)
	{
		m = me;
		v = va;
		sd = std::sqrt(va);
		cp = (ui)0;
	}
			
	template<class T>
	stats<T>::stats(long long s)
	{
		m = (T)0;
		v = (T)0;
		sd = (T)0;
		cp = (ui)0;
		ran.seed(s);
	}
	
	template<class T>
	stats<T>::stats(const T me,const T va,long long s)
	{
		m = me;
		v = va;
		sd = std::sqrt(va);
		cp = (ui)0;
		ran.seed(s);
	}
	
	template<class T>
	stats<T>::stats(const stats<T> &st)
	{
		m = st.m;
		v = st.v;
		sd = st.sd;
		cp = st.cp;
		ran = st.ran;
	}
	
	template<class T>
	stats<T>::~stats(void)
	{
		/*Do nothing*/
	}
	
	template<class T>
	stats<T> & stats<T>::operator=(const stats<T> &st)
	{
		m = st.m;
		v = st.v;
		sd = st.sd;
		cp = st.cp;
		ran = st.ran;
		return *this;
	}
	
	template<class T>
	stats<T> & stats<T>::mean(const T me)
	{
		m = me;
		return *this;
	}
	
	template<class T>
	inline
	stats<T> & stats<T>::var(const T va)
	{
		v = va;
		sd = std::sqrt(va);	
		return *this;
	} 
	
	template<class T>
	inline
	stats<T> & stats<T>::count(const unsigned int co)
	{
		cp = co;
		return *this;
	}
	
	template<class T>
	inline
	stats<T> & stats<T>::stdd(const T st)
	{
		sd = st;
		v = st*st;	
		return *this;
	}
	
	template<class T>
	inline
	T & stats<T>::mean(void)
	{
		return m;
	}
	
	template<class T>
	inline
	T & stats<T>::var()
	{
		return v;
	} 
	
	template<class T>
	inline
	unsigned int & stats<T>::count()
	{
		return cp;
	}
	
	template<class T>
	inline
	T & stats<T>::stdd()
	{
		return sd;
	}
	
	/*Return our random generator object reference*/
	template<class T>
	inline
	random_generator<T> & stats<T>::random(void)
	{
		return ran;
	}
	
	template<class T>
	inline
	stats<T> & stats<T>::update(const T x_new)
	{
		if(cp == 0) 
		{
			m = x_new; 
			++cp;
			v = (T) 0;
			sd = (T) 0;
			//std::cout << m << "\n";
			return *this;
		} 
		else if(cp >= 1)
		{
			/*Update the mean*/
			m = ((T)1/(T)(cp+1))*(((T)cp)*m + x_new);
			/*Update the variance*/
			//v = ((T)(cp)*v + (x_new*x_new - m*m))/(T)(cp+1);
			v = ((T)1/(T)(cp))*((T)(cp-1)*v + (x_new-m)*(x_new-m));
			/*Update the std deviation*/
			sd = std::sqrt(v);
			++cp;
			return *this;
		}
		else
		{
			std::cout << "Runtime Error\nSomething went wrong\n"; exit(1);
		}
	}
	
	template<class T>
	inline
	stats<T> & stats<T>::reset(void)
	{
		cp = (unsigned int)0;
		m = (T)0;
		v = (T)0;
		sd = (T)0;
		return *this;
	}
	
	template<class T>
	inline
	T stats<T>::confidance_interval(const unsigned int per)
	{
		/*Simple funciton that takes care of confidance intervals*/
		T z = (T)0.;
		/*All the special z - table values we could use*/
		switch(per)
		{
			case 99: z = (T)2.58; break;
			case 98: z = (T)2.33; break;
			case 95: z = (T)1.96; break;
			case 90: z = (T)1.64; break;
			case 85: z = (T)1.44; break;
			case 80: z = (T)1.28; break;
			default: std::cout << "Runtime Error\nNot a correct/standard confidance level\n"; exit(1);break;
		}
		//std::cout << "Confidance Interval " << per << "%: " << m - z*sd/(T)std::sqrt(cp-1) << " < µ < " << m + z*sd/(T)std::sqrt(cp-1) << "\n";
		return ((T)2.*z*sd/(T)std::sqrt((T)cp-(T)0));
	}
	
	template<class T>
	T stats<T>::expected(T (*sim)(random_generator<T>&),const long long N,const unsigned int per)
	{
		if(N < 2) {std::cout << "Runtime Error\nRealizations needs to be greater then 2\n"; exit(1);} else {}
		/*Iterate over the number of realizations*/
		for(long long i = 0; i < N; ++i)
		{
			update(sim(random())); //Update the stats with the new value
		}
		/*Generate a confidance interval*/
		T cI = confidance_interval(per);
		std::cout << std::setprecision(8) << "Expected: " << m << " Realizations: " << N << " Confidance Interval " << per << "%: " << m-cI << " < µ < " << m+cI << " Interval Length: " << cI << " MSE: " << v/cp << " RMSE: " << sd/std::sqrt(cp) << "\n"; 
		return m;
	}
	
	template<class T>
	T stats<T>::expected(T (*sim)(random_generator<T>&),const T tol,const unsigned int per)
	{
		T cI = tol+.1;
		/*Iterate until cI < tol -> Also includes a max number of iteratiosn just incase it is infinite loop*/
		while((cI > tol || cp < 100) && cp != 100000000000)
		{
			/*Update the statistics*/
			update(sim(random()));
			/*Generate the length of our interval to see if we have a good tolerance*/
			cI = confidance_interval(per);
		}
		std::cout << "Expected: " << m << " Realizations: " << cp << " Confidance Interval " << per << "%: " << m-cI << " < µ < " << m+cI << " Interval Length: " << cI << " MSE: " << v/cp << " RMSE: " << sd/std::sqrt(cp) << "\n"; 
		return m;
	}
	
	template<class T>
	T stats<T>::var(T (*sim)(random_generator<T>&),const long long N,const long long boot,const int per)
	{
		/*Save our realizations*/
		matrix<T> data = matrix<T>::zeros((ui)N,1);
		/*Save our Bootstap estimates*/
		matrix<T> Y = matrix<T>::zeros((ui)boot,1);
		/*Run N realizations*/
		for(long long i = 0; i < N; ++i)
		{
			/*Generate a realization*/
			T X = sim(random());
			/*Save the data from the realization*/
			data((ui)i) = X;
			/*Update the other results*/
			update(X);
		}
		/*Save our variance estimate*/
		const T cp_var = v;
		
		/*Reset results to use to estimate the estimator error*/
		reset();
		
		/*Now to update the MSE of the Variance estimator*/
		for(int i = 0; i < (int)boot; ++i)
		{
			/*Random sampling n times for our estimator*/
			for(int j = 0; j < (int)N; ++j)
			{
				/*Pick random sample from our data vector and update our stat values*/
				update(data((unsigned int)random().uniform((T)0.,(T)(N))));
			}
			/*Update our MSE estimate*/
			Y(i) = (v - cp_var)*(v - cp_var);
			reset();
		}
		
		/*Calculate the mean (Do it here to avoid name conflicts)*/
		T mean_Y = (T)0;
		for(int i = 0; i < boot; ++i)
		{
			mean_Y += Y((ui)i);
		}
		/*MSE*/
		mean_Y /= boot;
		
		/*Get the confidance level of our estimator using MSE*/
		T z = (T)0.;
		/*All the special z - table values we could use*/
		switch(per)
		{
			case 99: z = (T)2.58; break;
			case 98: z = (T)2.33; break;
			case 95: z = (T)1.96; break;
			case 90: z = (T)1.64; break;
			case 85: z = (T)1.44; break;
			case 80: z = (T)1.28; break;
			default: std::cout << "Runtime Error\nNot a correct/standard confidance level\n"; exit(1);break;
		}
		
		/*Return results*/
		std::cout << "Var: " << cp_var << " Number of realizations: " << N << " Confidance level " << per << "%: " << cp_var - z*std::sqrt(mean_Y) <<  " < σ^2 < " << cp_var + z*std::sqrt(mean_Y) << " Length of interval: " << 2.*z*std::sqrt(mean_Y) << " MSE: " << mean_Y  << " RMSE: "<< std::sqrt(mean_Y) << "\n";
		return cp_var;
	}
	
	template<class T>
	T stats<T>::expected_bootstrap(T (*sim)(random_generator<T>&),const long long N,const long long boot,const int per)
	{
		/*Save our realizations*/
		matrix<T> data = matrix<T>::zeros((ui)N,1);
		/*Save our Bootstap estimates*/
		matrix<T> Y = matrix<T>::zeros((ui)boot,1);
		/*Run N realizations*/
		for(long long i = 0; i < N; ++i)
		{
			/*Generate a realization*/
			T X = sim(random());
			/*Save the data from the realization*/
			data((ui)i) = X;
			/*Update the other results*/
			update(X);
		}
		/*Save our variance estimate*/
		const T cp_mean = m;
		
		/*Reset results to use to estimate the estimator error*/
		reset();
		
		/*Now to update the MSE of the Variance estimator*/
		for(int i = 0; i < (int)boot; ++i)
		{
			/*Random sampling n times for our estimator*/
			for(int j = 0; j < (int)N; ++j)
			{
				/*Pick random sample from our data vector and update our stat values*/
				update(data((unsigned int)random().uniform((T)0.,(T)(N))));
			}
			/*Update our MSE estimate*/
			Y(i) = (m - cp_mean)*(m - cp_mean);
			reset();
		}
		
		/*Calculate the mean (Do it here to avoid name conflicts)*/
		T mean_Y = (T)0;
		for(int i = 0; i < boot; ++i)
		{
			mean_Y += Y((ui)i);
		}
		/*MSE*/
		mean_Y /= boot;
		
		/*Get the confidance level of our estimator using MSE*/
		T z = (T)0.;
		/*All the special z - table values we could use*/
		switch(per)
		{
			case 99: z = (T)2.58; break;
			case 98: z = (T)2.33; break;
			case 95: z = (T)1.96; break;
			case 90: z = (T)1.64; break;
			case 85: z = (T)1.44; break;
			case 80: z = (T)1.28; break;
			default: std::cout << "Runtime Error\nNot a correct/standard confidance level\n"; exit(1);break;
		}
		
		/*Return results*/
		std::cout << "Expected: " << cp_mean << " Number of realizations: " << N << " Confidance level " << per << "%: " << cp_mean - z*std::sqrt(mean_Y) <<  " < µ < " << cp_mean + z*std::sqrt(mean_Y) << " Length of interval: " << 2.*z*std::sqrt(mean_Y) << " MSE: " << mean_Y  << " RMSE: "<< std::sqrt(mean_Y) << "\n";
		return cp_mean;
	}
	
	template<class T>
	T stats<T>::greater(T (*sim)(random_generator<T>&),const long long N,const unsigned int per,const T c)
	{
		if(N < 2) {std::cout << "Runtime Error\nRealizations needs to be greater then 2\n"; exit(1);} else {}
		/*Create our new indicator function*/
		auto indicator = [](T (*sim)(random_generator<T>&),random_generator<T> &r,T c)
		{
			T x = (T)sim(r);
			if(x > c) {return (T)1.;}
			else {return (T)0;}
		};
		
		/*Iterate over the number of realizations*/
		for(long long i = 0; i < N; ++i)
		{
			update(indicator(sim,random(),c)); //Update the stats with the new value
		}
		/*Generate a confidance interval*/
		T cI = confidance_interval(per);
		std::cout << "P{X > c}: " << m << " Realizations: " << N << " Confidance Interval " << per << "%: " << m-cI << " < P{X > c} < " << m+cI << " Interval Length: " << cI << " MSE: " << v/cp << " RMSE: " << sd/std::sqrt(cp) << "\n"; 
		return m;
	}
	
	template<class T>
	T stats<T>::less(T (*sim)(random_generator<T>&),const long long N,const unsigned int per,const T c)
	{
		if(N < 2) {std::cout << "Runtime Error\nRealizations needs to be greater then 2\n"; exit(1);} else {}
		/*Create our new indicator function*/
		auto indicator = [](T (*sim)(random_generator<T>&),random_generator<T> &r,T c)
		{
			T x = (T)sim(r);
			if(x < c) {return (T)1.;}
			else {return (T)0;}
		};
		
		/*Iterate over the number of realizations*/
		for(long long i = 0; i < N; ++i)
		{
			update(indicator(sim,random(),c)); //Update the stats with the new value
		}
		/*Generate a confidance interval*/
		T cI = confidance_interval(per);
		std::cout << "P{X < c}: " << m << " Realizations: " << N << " Confidance Interval " << per << "%: " << m-cI << " < P{X < c} < " << m+cI << " Interval Length: " << cI << " MSE: " << v/cp << " RMSE: " << sd/std::sqrt(cp) << "\n"; 
		return m;
	}
	
	template<class T>
	T stats<T>::greater_equal(T (*sim)(random_generator<T>&),const long long N,const unsigned int per,const T c)
	{
		if(N < 2) {std::cout << "Runtime Error\nRealizations needs to be greater then 2\n"; exit(1);} else {}
		/*Create our new indicator function*/
		auto indicator = [](T (*sim)(random_generator<T>&),random_generator<T> &r,T c)
		{
			T x = (T)sim(r);
			if(x >= c) {return (T)1.;}
			else {return (T)0;}
		};
		
		/*Iterate over the number of realizations*/
		for(long long i = 0; i < N; ++i)
		{
			update(indicator(sim,random(),c)); //Update the stats with the new value
		}
		/*Generate a confidance interval*/
		T cI = confidance_interval(per);
		std::cout << "P{X >= c}: " << m << " Realizations: " << N << " Confidance Interval " << per << "%: " << m-cI << " < P{X >= c} < " << m+cI << " Interval Length: " << cI << " MSE: " << v/cp << " RMSE: " << sd/std::sqrt(cp) << "\n"; 
		return m;
	}
	
	template<class T>
	T stats<T>::less_equal(T (*sim)(random_generator<T>&),const long long N,const unsigned int per,const T c)
	{
		if(N < 2) {std::cout << "Runtime Error\nRealizations needs to be greater then 2\n"; exit(1);} else {}
		/*Create our new indicator function*/
		auto indicator = [](T (*sim)(random_generator<T>&),random_generator<T> &r,T c)
		{
			T x = (T)sim(r);
			if(x <= c) {return (T)1.;}
			else {return (T)0;}
		};
		
		/*Iterate over the number of realizations*/
		for(long long i = 0; i < N; ++i)
		{
			update(indicator(sim,random(),c)); //Update the stats with the new value
		}
		/*Generate a confidance interval*/
		T cI = confidance_interval(per);
		std::cout << "P{X <= c}: " << m << " Realizations: " << N << " Confidance Interval " << per << "%: " << m-cI << " < P{X <= c} < " << m+cI << " Interval Length: " << cI << " MSE: " << v/cp << " RMSE: " << sd/std::sqrt(cp) << "\n"; 
		return m;
	}
	
	template<class T>
	T stats<T>::equal(T (*sim)(random_generator<T>&),const long long N,const unsigned int per,const T c)
	{
		if(N < 2) {std::cout << "Runtime Error\nRealizations needs to be greater then 2\n"; exit(1);} else {}
		/*Create our new indicator function*/
		auto indicator = [](T (*sim)(random_generator<T>&),random_generator<T> &r,T c)
		{
			T x = (T)sim(r);
			if(x == c) {return (T)1.;}
			else {return (T)0;}
		};
		
		/*Iterate over the number of realizations*/
		for(long long i = 0; i < N; ++i)
		{
			update(indicator(sim,random(),c)); //Update the stats with the new value
		}
		/*Generate a confidance interval*/
		T cI = confidance_interval(per);
		std::cout << "P{X == c}: " << m << " Realizations: " << N << " Confidance Interval " << per << "%: " << m-cI << " < P{X == c} < " << m+cI << " Interval Length: " << cI << " MSE: " << v/cp << " RMSE: " << sd/std::sqrt(cp) << "\n"; 
		return m;
	}
	
	template<class T>
	T stats<T>::expected(T (*sim)(random_generator<T>&),const long long N,const unsigned int per,T (*sim2)(random_generator<T>&),const T mu_y)
	{
		/*Vectors to store the first 100 values for both simulations*/
		matrix<T> X = matrix<T>::zeros(100,1);
		matrix<T> Y = matrix<T>::zeros(100,1);
		
		/*Save the mean*/
		T x_mean = (T)0.,y_mean = (T)0.;
		/*Start the simulation*/
		for(int i = 0; i < X.numel(); ++i)
		{
			/*Save the seed*/
			long long seed = ran.seed();
			/*Run the simulation*/
			X(i) = sim(ran);
			/*Reset the seed*/
			ran.seed(seed);
			Y(i) = sim2(ran);
			//std::cout << Y(i) << "\n";
			
			/*Get the mean*/
			x_mean += X(i);
			y_mean += Y(i);
		}
		/*Get the average value*/
		x_mean /= (T)X.numel();
		y_mean /= (T)Y.numel();
		
		/*Calculate the covariance and variance*/
		T y_var = (T)0.,cv = (T)0.;
		for(int i = 0; i < X.numel(); ++i)
		{
			y_var += (Y(i) - y_mean)*(Y(i) - y_mean);
			cv += (Y(i) - y_mean)*(X(i) - x_mean);
		}
		
		/*Update what the most optimal c is*/
		T c = -cv/y_var;
		//std::cout << "\n" << cv/(T)(100-1) << "\n" << y_var/(T)(100-1) << "\n" << c << "\n";
		
		/*Update the previous results*/
		for(int i = 0; i < X.numel(); ++i)
		{
			update(X(i) + c*(Y(i) - mu_y));
		}
		
		/*Emply our new linear combination to reduce the variance*/
		for(long long i = (ll)X.numel(); i < N; ++i)
		{
			/*Save the seed*/
			long long seed = ran.seed();
			/*Update the result we know*/
			T Z = sim(ran);
			/*Reset the seed to keep the corralation*/
			ran.seed(seed);
			/*Update the new Random variable*/
		    Z += c*(sim2(ran) - mu_y);
			update(Z);
		}
		/*Generate a confidance interval*/
		T cI = confidance_interval(per);
		std::cout << "Expected: " << m << " Realizations: " << N << " Confidance Interval " << per << "%: " << m-cI << " < µ < " << m+cI << " Interval Length: " << cI << " MSE: " << v/cp << " RMSE: " << sd/std::sqrt(cp) << "\n";
		std::cout << "C_opt: " << c << "\n";
		return m;
	}
	
	template<class T>
	T stats<T>::expected(T (*sim)(random_generator<T>&),const long long N,const unsigned int per,T (*sim2)(random_generator<T>&),const T mu_y,T (*sim3)(random_generator<T>&),const T mu_r)
	{
		/*Vectors to store the first 100 values for the 3 simulations*/
		matrix<T> X = matrix<T>::zeros(100,1);
		matrix<T> Y = matrix<T>::zeros(100,1);
		matrix<T> R = matrix<T>::zeros(100,1);
		
		/*Store the covariance matrix and b vector for our optimal cs*/
		matrix<T> C = matrix<T>::zeros(2,2);
		matrix<T> b = matrix<T>::zeros(2,1);
		
		/*Mean estimator for the simulations*/
		T mean_x = (T)0, mean_y = (T)0, mean_r = (T)0;
		for(int i = 0; i < X.numel(); ++i)
		{
			/*Save the seed to reuse with corralation (This process can be optimized if I had a random generator that didnt update)*/
			long long seed = ran.seed();
			/*Update the simulated values*/
			X(i) = sim(ran);
			ran.seed(seed);
			Y(i) = sim2(ran);
			ran.seed(seed);
			R(i) = sim3(ran);
			
			/*Update the mean estimate*/
			mean_x += X(i);
			mean_y += Y(i);
			mean_r += R(i);
		}
		/*Update the mean values*/
		mean_x /= (T)X.numel();
		mean_y /= (T)X.numel();
		mean_r /= (T)X.numel();
		
		/*Generate covariances and variances*/
		for(int i = 0; i < X.numel(); ++i)
		{
			C(0,0) += (Y(i) - mean_y)*(Y(i) - mean_y);
			C(1,1) += (R(i) - mean_r)*(R(i) - mean_r);
			C(0,1) += (Y(i) - mean_y)*(R(i) - mean_r);
			C(1,0) = C(0,1);
			b(0) += -(X(i) - mean_x)*(Y(i) - mean_y);
			b(1) += -(X(i) - mean_x)*(R(i) - mean_r);
		}
		/*Solve for optimal c*/
		matrix<T> c = gauss<T>::solve(C,b);
		
		/*Redo the previous results with our new optimal C*/
		for(int i = 0; i < X.numel(); ++i)
		{
			update(X(i) + c(0)*(Y(i) - mu_y) + c(1)*(R(i) - mu_r));
		}
		
		/*Redo the remaining realizations*/
		for(long long i = (ll)X.numel(); i < N; ++i)
		{
			/*Save the seed so the results will have the same corralation structure*/
			long long seed = ran.seed();
			T x = sim(ran);
			ran.seed(seed);
			T y = sim2(ran) - mu_y;
			ran.seed(seed);
			T r = sim3(ran) - mu_r;
			
			/*Put together and update our realization*/
			update(x + c(0)*y + c(1)*r);
		}
		/*Generate a confidance interval*/
		T cI = confidance_interval(per);
		std::cout << std::setprecision(8) << "Expected: " << m << " Realizations: " << N << " Confidence Interval " << per << "%: " << m-cI << " < µ < " << m+cI << " Interval Length: " << cI << " MSE: " << v/cp << " RMSE: " << sd/std::sqrt(cp) << "\n"; 
		std::cout << "C_opt: " << "(" << c(0) << "," << c(1) << ")\n";
		return m;
	}
	
	/*The expected value with strata uniforms*/
	template<class T>
	T stats<T>::expected_strata(T (*sim)(random_generator<T>&,T,T),const long long N,const unsigned int per,const unsigned int k)
	{
		if(N < 2) {std::cout << "Runtime Error\nRealizations needs to be greater then 2\n"; exit(1);} else {}
		/*Iterate over the number of realizations and get one stata per realization (This ignores the end points)*/
		if(N%k != 0) {std::cout << "Runtime Error\nk and N must have no remainder since it isnt programed\n"; exit(1);}
		
		/*Create a new stats class to store the other means*/
		stats<T> st; st.reset();
		
		/*Number of realizations each strata*/
		unsigned int n_i = (unsigned int) N/k;
		
		/*Save the variance of each strata*/
		T var_s = (T)0;
 		
		/*Iterate over all the strata*/
		for(int j = 0; j < k; ++j)
		{
			/*Do realizations in strata*/
			for(int n = 0; n < n_i; ++n)
			{
				st.update(sim(random(),(T)j/(T)k,(T)(j+1)/(T)k));
			}
			update(st.mean()); var_s += st.var();
			st.reset();
		}
		/*Update our variance for our estimator*/
		v = var_s/((T)k);
	    sd = std::sqrt(v);
		
		/*Generate a confidance interval*/
		T cI = confidance_interval(per);
		std::cout << std::setprecision(8) << "Expected: " << m << " Realizations: " << N << " Confidence Interval " << per << "%: " << m-cI << " < µ < " << m+cI << " Interval Length: " << cI << " MSE: " << v/(T)cp << " RMSE: " << sd/(T)std::sqrt(cp) << "\n"; 
		return m;
	}
	
	template<class T>
	T stats<T>::gaussian(const T x)
	{
		return ((T)1/(std::sqrt(2*pi)))*std::exp(-(x*x)/((T)2));
	}
	
	template<class T>
	T stats<T>::gaussian(const T mu,const T sigma,const T x)
	{
		return ((T)1/(std::sqrt(2*pi*sigma)))*std::exp(-((x-mu)*(x-mu))/((T)2*sigma));
	}
	
	template<class T>
	T stats<T>::exponential(const T lambda,const T x)
	{
		return lambda*std::exp(-lambda*x);
	}
	
	template<class T>
	T stats<T>::uniform(void)
	{
		return (T)1;
	}
	
	template<class T>
	T stats<T>::uniform(T a, T b)
	{
		return (T)1/(b - a);
	}
	
}