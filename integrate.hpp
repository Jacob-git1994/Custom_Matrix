/*
	Jacob Dresher
	Trapazoidal Numerical integration using Richardson 
	extrapolation to increase the order of the method many 
	levels to increase numerical appoximation

	Code Written using Object oriented programing and allows customizations also with fopenmp to speed up large numerical grids
	Note: If the grid or numbers are too large you will get a meaningless result
	
	Later versions will incorporate mulidimentional grids with equal spacing

	Note: Periodic functions spectrally converge
*/

#pragma once
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <chrono>

/*Some useful constant for example problems*/
#define Pi 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384460955058223172535940812848111745028410270

namespace nicole
{
	template<class T>
	class integrate
	{
	private:
		/*void one_step_trap_1D(type func(type),lower bound,upper bound,gridspacing)*/
		T trap_1D(T (*)(const T),const T,const T,const unsigned int);
		/*2D trap*/
		T trap_2D(T (*)(const T,const T),const T,const T,const T,const T,const unsigned int);
		/*Number of dimentions to integrate over*/
		unsigned int DIM;
		/*Specify the inital ∆x*/
		unsigned int inital_step;
		/*Number of levels we want to extrapolate the solution to*/
		unsigned int num_levels;
		/*Error tolerence between the appoximations*/
		T tol;
		/*Generate a lower triangular matrix to store the appoximations*/
		T **Lower;
		/*Generate an array based of the dimention bounds*/
		T *bounds;
	public:
		/*Void integrate(void)*/
		integrate(void);
		/*1D Void integrate(lower,upper,steps,tolerence,number of levels)*/
		integrate(const T,const T,const unsigned int,const T,const unsigned int);
		/*Integrate without worry about the number of levels:Note leads to more inaccuracy*/
		integrate(const T,const T,const unsigned int);
		
		/*2D Void integrate*/
		integrate(const T,const T,const T,const T,const unsigned int,const T,const unsigned int);
		/*2D integrate without worring about the number of levels:Note leads to more inaccuracy*/
		integrate(const T,const T,const T,const T,const unsigned int);
		
		/*Destructor*/
		~integrate(void);
		/*Higher dimentional implimentations coming later grid spacing must be the same*/
		/*3D Void integrate(lower1,upper1,lower2,upper2,lower3,upper3,step size x,step size y,step size z,tolerence,number of levels)*/
		
		/*Functions for 1D cases*/
		integrate<T> & set_bounds(const T,const T);
		integrate<T> & set_bounds(const T,const T,const unsigned int);
		
		/*Functions for 2D case*/
		integrate<T> & set_bounds(const T,const T,const T,const T);
		integrate<T> & set_bounds(const T,const T,const T,const T,const unsigned int);
		
		/*Set the number of appoximation levels*/
		integrate<T> & set_num_grid_points(const unsigned int);
		integrate<T> & set_num_levels(const unsigned int);
		integrate<T> & set_tol(const T);
		
		/*Prints the Lower Triangular Matrix*/
		integrate<T> & print_lower(void);
		integrate<T> & print_dim(void) {std::cout << "Dimentions: " << DIM << "\n"; return *this;}
		
		/*Solve*/
		T solve(T (*)(const T));
		T solve(T (*)(const T,const T));
	};
	
	template<class T>
	integrate<T>::integrate(void)
	{
		DIM = (unsigned int)0;
		inital_step = (unsigned int)0;
		num_levels = (unsigned int)0;
		tol = (T).0001;
		
		Lower = nullptr;
		bounds = nullptr;
	}
	
	template<class T>
	integrate<T>::integrate(const T L,const T U,const unsigned int step,const T TOL,const unsigned int Lev)
	{
		/*Initalize the levels*/
		DIM = (unsigned int)1;
		if(TOL <= 0) {std::cout<<"Runtime Error\nTol cant be ≤ 0\nExiting\n"; exit(1);} else {tol = TOL;}
		if(Lev <= 0) {std::cout<<"Runtime Error\nLevels cant be ≤ 0\nExiting\n"; exit(1);} else {num_levels = Lev;}
		if(step <= 0) {std::cout<<"Runtime Error\nstep size cant be ≤ 0\nExiting\n"; exit(1);} else {inital_step = step;}
		
		/*Initalize the Lower Triangular matrix*/
		Lower = new T*[num_levels]; if(!Lower) {std::cout << "Runtime Error\nMemory Allocation Failure\nExiting\n";exit(1);} else {}
		for (int i = 0; i < num_levels; i++)
		{
			Lower[i] = new T[i+1]; if(!Lower[i]) {std::cout << "Runtime Error\nMemory Allocation Failure\nExiting\n";exit(1);} else {}
			for(int j = 0; j <= i; j++)
			{
				Lower[i][j] = (T)0;
			}
		}
		
		/*Initalize the Array for bounds*/
		/*This integration constructor only supports 1D*/
		bounds = new T[2]; if(!bounds) {std::cout << "Runtime Error\nMemory Allocation Failure\nExiting\n";exit(1);} else {}
		bounds[0] = L;
		bounds[1] = U;
		//std::cout << bounds[0] << " " << bounds[1] << "\n" << inital_step_size << "\n" << num_levels << "\n";
	}
	
	template<class T>
	integrate<T>::integrate(const T L,const T U,const unsigned int step)
	{
		num_levels = (unsigned int)1;
		tol = .001; /*Does not matter here*/
		if(step <= 0) {std::cout<<"Runtime Error\nstep size cant be ≤ 0\nExiting\n"; exit(1);} else {inital_step = step;}
		
		/*Initalize the Lower Triangular matrix*/
		Lower = new T*[num_levels]; if(!Lower) {std::cout << "Runtime Error\nMemory Allocation Failure\nExiting\n";exit(1);} else {}
		for (int i = 0; i < num_levels; i++)
		{
			Lower[i] = new T[i+1]; if(!Lower[i]) {std::cout << "Runtime Error\nMemory Allocation Failure\nExiting\n";exit(1);} else {}
			for(int j = 0; j <= i; j++)
			{
				Lower[i][j] = (T)0;
			}
		}
		/*Initalize the Array for bounds*/
		/*This integration constructor only supports 1D*/
		bounds = new T[2]; if(!bounds) {std::cout << "Runtime Error\nMemory Allocation Failure\nExiting\n";exit(1);} else {}
		bounds[0] = L;
		bounds[1] = U;
	}
	
	/*2D setting up of our problem*/
	template<class T>
	integrate<T>::integrate(const T L1,const T U1,const T L2,const T U2,const unsigned int step,const T TOL,const unsigned int Lev)
	{
		/*Initalize the levels*/
		DIM = (unsigned int)2;
		if(TOL <= 0) {std::cout<<"Runtime Error\nTol cant be ≤ 0\nExiting\n"; exit(1);} else {tol = TOL;}
		if(Lev <= 0) {std::cout<<"Runtime Error\nLevels cant be ≤ 0\nExiting\n"; exit(1);} else {num_levels = Lev;}
		if(step <= 0) {std::cout<<"Runtime Error\nstep size cant be ≤ 0\nExiting\n"; exit(1);} else {inital_step = step;}
		
		/*Initalize the Lower Triangular matrix*/
		Lower = new T*[num_levels]; if(!Lower) {std::cout << "Runtime Error\nMemory Allocation Failure\nExiting\n";exit(1);} else {}
		for (int i = 0; i < num_levels; i++)
		{
			Lower[i] = new T[i+1]; if(!Lower[i]) {std::cout << "Runtime Error\nMemory Allocation Failure\nExiting\n";exit(1);} else {}
			for(int j = 0; j <= i; j++)
			{
				Lower[i][j] = (T)0;
			}
		}
		
		/*Initalize the Array for bounds*/
		/*This integration constructor only supports 1D*/
		bounds = new T[4]; if(!bounds) {std::cout << "Runtime Error\nMemory Allocation Failure\nExiting\n";exit(1);} else {}
		bounds[0] = L1;
		bounds[1] = U1;
		bounds[2] = L2;
		bounds[3] = U2;
		//std::cout << bounds[0] << " " << bounds[1] << "\n" << inital_step_size << "\n" << num_levels << "\n";
	}

	template<class T>
	integrate<T>::~integrate(void)
	{
		delete[] bounds;
		for(int i = num_levels-1; i >= 0; i--)
		{
			delete[] Lower[i];
		}
		delete[] Lower;
		std::cout << "Finished Integration!\n";
	}

	/*Functions to set different condtions in 1D*/
	template<class T>
	integrate<T> & integrate<T>::set_bounds(const T L,const T U)
	{
		if(bounds == nullptr) {/*Do nothing*/} else {delete[] bounds;}
		bounds = new T[2]; if(!bounds) {std::cout << "Runtime Error\nMemory Allocation Failure\nExiting\n";exit(1);} else {}
		bounds[0] = L;
		bounds[1] = U;
		return *this;
	}

	template<class T>
	integrate<T> & integrate<T>::set_bounds(const T L,const T U,const unsigned int step)
	{
		if(step <= 0) {std::cout<<"Runtime Error\nstep size cant be ≤ 0\nExiting\n"; exit(1);} else {inital_step = step;}
		if(bounds == nullptr) {/*Do nothing*/} else {delete[] bounds;}
		bounds = new T[2]; if(!bounds) {std::cout << "Runtime Error\nMemory Allocation Failure\nExiting\n";exit(1);} else {}
		bounds[0] = L;
		bounds[1] = U;
		return *this;
	}
	
	template<class T>
	integrate<T> & integrate<T>::set_bounds(const T L1,const T U1,const T L2,const T U2)
	{
		
		if(bounds == nullptr) {/*Do nothing*/} else {delete[] bounds;}
		bounds = new T[4]; if(!bounds) {std::cout << "Runtime Error\nMemory Allocation Failure\nExiting\n";exit(1);} else {}
		bounds[0] = L1;
		bounds[1] = U1;
		bounds[2] = L2;
		bounds[3] = U2;
		return *this;
	}
	
	template<class T>
	integrate<T> & integrate<T>::set_bounds(const T L1,const T U1,const T L2,const T U2,const unsigned int step)
	{
		if(step <= 0) {std::cout<<"Runtime Error\nstep size cant be ≤ 0\nExiting\n"; exit(1);} else {inital_step = step;}
		if(bounds == nullptr) {/*Do nothing*/} else {delete[] bounds;}
		bounds = new T[4]; if(!bounds) {std::cout << "Runtime Error\nMemory Allocation Failure\nExiting\n";exit(1);} else {}
		bounds[0] = L1;
		bounds[1] = U1;
		bounds[2] = L2;
		bounds[3] = U2;
		return *this;
	}

	template<class T>
	integrate<T> & integrate<T>::set_num_grid_points(const unsigned int step)
	{
		if(step <= 0) {std::cout<<"Runtime Error\nstep size cant be ≤ 0\nExiting\n"; exit(1);} else {inital_step = step;}
		return *this;
	}

	template<class T>
	integrate<T> & integrate<T>::set_num_levels(const unsigned int Lev)
	{
		/*Delete the inital object data*/
		for(int i = num_levels-1; i >= 0; i--)
		{
			delete[] Lower[i];
		}
		delete[] Lower;
		
		/*Update the size*/
		if(Lev <= 0) {std::cout<<"Runtime Error\nLevels cant be ≤ 0\nExiting\n"; exit(1);} else {num_levels = Lev;}
		
		/*Initalize the Lower Triangular matrix again*/
		Lower = new T*[num_levels]; if(!Lower) {std::cout << "Runtime Error\nMemory Allocation Failure\nExiting\n";exit(1);} else {}
		for (int i = 0; i < num_levels; i++)
		{
			Lower[i] = new T[i+1]; if(!Lower[i]) {std::cout << "Runtime Error\nMemory Allocation Failure\nExiting\n";exit(1);} else {}
			for(int j = 0; j <= i; j++)
			{
				Lower[i][j] = (T)0;
			}
		}
		return *this;
	}

	template<class T>
	integrate<T> & integrate<T>::set_tol(const T TOL)
	{
		if(TOL <= 0) {std::cout<<"Runtime Error\nTol cant be ≤ 0\nExiting\n"; exit(1);} else {tol = TOL;}
		return *this;
	}
	
	/*Method to analize the lower traiangular matrix if the method did not converge to the right tolerance*/
	template<class T>
	integrate<T> & integrate<T>::print_lower(void)
	{
		for(int i = 0; i < num_levels; i++)
		{
			for(int j = 0; j <= i; j++)
			{
				std::cout << Lower[i][j] << "\t";
			}
			std::cout << "\n";
		}
		return *this;
	}

	/*Impliment the trapazoid rule*/
	template<class T>
	T integrate<T>::trap_1D(T (*f)(const T),const T L,const T U,const unsigned int step)
	{
		/*Set up grid spacing*/
		const T h = ((T)1)/((T)step);
		/*Initalize the sum*/
		T sum = (T).5*(f((U-L)*(T)0 + L) + f((U-L)*(T)1) + L);
		
		/*Iterate over the internal points maped from x -> f(x)*/
#pragma omp parallel for reduction(+:sum) schedule(static) 
		for(int i = 1; i < step; i++)
		{
			sum += f((U-L)*h*(T)i + L);
		}
		/*Multiply the ∆x to the sum*/
		sum *= h;
		/*Return the sum*/
		return (U-L)*sum;
	}

	template<class T>
	T integrate<T>::solve(T (*f)(const T))
	{
		auto t1 = std::chrono::high_resolution_clock::now();
		//T fin = trap_1D(f,bounds[0],bounds[1],inital_step);
		
		/*Compute the inital step*/
		Lower[0][0] = trap_1D(f,bounds[0],bounds[1],inital_step);
	
		for(int i = 1; i < num_levels; i++)
		{
			Lower[i][0] = trap_1D(f,bounds[0],bounds[1],std::pow((T)2,(T)i)*inital_step);
			for(int j = 1; j <= i; j++)
			{
				Lower[i][j] = ((std::pow((T)4,(T)j))*Lower[i][j-1] - Lower[i-1][j-1])/(std::pow((T)4,(T)j) - (T)1);
			}
			
			/*Check the Error Tol on the diagonal*/
			if(std::fabs(Lower[i][i] - Lower[i-1][i-1])/std::fabs(Lower[i][i]) <= tol) 
			{
				std::cout << "Solved in: " << std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - t1).count() << " seconds\n" \
					<< "Reletive Error: " << std::fabs(Lower[i][i] - Lower[i-1][i-1])/std::fabs(Lower[i][i]) << "\n";
				return Lower[i][i];
			}
			else {} 
		}
		
		std::cout << "Solved in " << std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - t1).count() << " seconds\n";
		if(num_levels > 1)
		{
			std::cout << "We werent able to find a solution in the desirered tolerance\nReturning best estimate:\n";
			std::cout << "Reletive Error: " << std::fabs(Lower[num_levels-1][num_levels-1] - Lower[num_levels-1][num_levels-2])/std::fabs(Lower[num_levels-1][num_levels-1]) << "\n";
			
		}
		else 
		{
			/*Do nothing*/
		}
		return Lower[num_levels-1][num_levels-1];
	}
	
	/*2D intergration*/
	template<class T>
	T integrate<T>::trap_2D(T (*f)(const T,const T),const T L1,const T U1,const T L2,const T U2,const unsigned int step)
	{
		/*Set up grid spacing*/
		const T h = ((T)1)/((T)step);
		/*Initalize where we save the sum*/
		T sum = (T)0;
		
		/*Sum the entire boundary >> No need to parallize this loop*/
		for(int i = 0; i <= step; i++)
		{
			sum += f(L1,(U2-L2)*h*i + L2) + f(U1,(U2-L2)*h*i + L2) + f((U1-L1)*h*i + L1,L2) + f((U1-L1)*h*i + L1,U2);
		}
		/*Multiply by .5*/
		sum *= (T).5;
		
		/*Sum the internal grid*/
#pragma omp parallel for reduction(+:sum) schedule(static)
		for(int i = 1; i < step; i++)
		{
			for(int j = 1; j < step; j++)
			{
				sum += f((U1-L1)*h*(T)i + L1,(U2-L2)*h*(T)j + L2);
			}
		}
		
		/*Multiply by the distances*/
		sum *= h*h;
		return (U1-L1)*(U2-L2)*sum;
	}
	
	/*Integrate over 2D bounds*/
	template<class T>
	T integrate<T>::solve(T (*f)(const T,const T))
	{
		auto t1 = std::chrono::high_resolution_clock::now();
		//T fin = trap_1D(f,bounds[0],bounds[1],inital_step);
		
		/*Compute the inital step*/
		Lower[0][0] = trap_2D(f,bounds[0],bounds[1],bounds[2],bounds[3],inital_step);
	
		for(int i = 1; i < num_levels; i++)
		{
			Lower[i][0] = trap_2D(f,bounds[0],bounds[1],bounds[2],bounds[3],std::pow((T)2,(T)i)*inital_step);
			for(int j = 1; j <= i; j++)
			{
				Lower[i][j] = ((std::pow((T)4,(T)j))*Lower[i][j-1] - Lower[i-1][j-1])/(std::pow((T)4,(T)j) - (T)1);
			}
			
			/*Check the Error Tol on the diagonal*/
			if(std::fabs(Lower[i][i] - Lower[i-1][i-1])/std::fabs(Lower[i][i]) <= tol) 
			{
				std::cout << "Solved in: " << std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - t1).count() << " seconds\n" \
					<< "Reletive Error: " << std::fabs(Lower[i][i] - Lower[i-1][i-1])/std::fabs(Lower[i][i]) << "\n";
				return Lower[i][i];
			}
			else {} 
		}
		
		std::cout << "Solved in " << std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - t1).count() << " seconds\n";
		if(num_levels > 1)
		{
			std::cout << "We werent able to find a solution in the desirered tolerance\nReturning best estimate:\n";
			std::cout << "Reletive Error: " << std::fabs(Lower[num_levels-1][num_levels-1] - Lower[num_levels-2][num_levels-2])/std::fabs(Lower[num_levels-1][num_levels-1]) << "\n";
			
		}
		else 
		{
			/*Do nothing*/
		}
		return Lower[num_levels-1][num_levels-1];
	}	
}