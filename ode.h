#pragma once
#include "matrix.h"
#include "gauss.h"
#include <stdlib.h>
#include <iostream>
#include <stdio.h>

namespace nicole
{
	/*When declaring the function vector make sure to create a matrix solution and a copy to allow the methods to work! Also make sure it is a const reference!*/
	template<class T>
	class ode
	{
	public:
		/*Forward Euler solvers ->> not as accurate but easy to impliment*/
		static matrix<T> euler(const matrix<T>&,matrix<T> (*)(const matrix<T>&,T),T,T,unsigned int);
		static matrix<T> euler(const matrix<T>&,matrix<T> (*)(const matrix<T>&),T,T,unsigned int);
		static matrix<T> euler(const matrix<T>&,matrix<T> (*)(const matrix<T>&,T),T,T,char [],unsigned int);
		static matrix<T> euler(const matrix<T>&,matrix<T> (*)(const matrix<T>&),T,T,char [],unsigned int);
		/*Forward euler with banded matrices that dont depend on time*/
		static matrix<T> euler(const matrix<T>&,const matrix<T>&,matrix<T> (*)(T),T,T,unsigned int,unsigned int); 
		static matrix<T> euler(const matrix<T>&,const matrix<T>&,matrix<T> (*)(T),T,T,char [],unsigned int,unsigned int);
		
		/*Runge_kutta second order or midpoint*/
		static matrix<T> rk2(const matrix<T>&,matrix<T> (*)(const matrix<T>&,T),T,T,unsigned int);
		static matrix<T> rk2(const matrix<T>&,matrix<T> (*)(const matrix<T>&),T,T,unsigned int);
		static matrix<T> rk2(const matrix<T>&,matrix<T> (*)(const matrix<T>&,T),T,T,char [],unsigned int);
		static matrix<T> rk2(const matrix<T>&,matrix<T> (*)(const matrix<T>&),T,T,char [],unsigned int);
		
		/*Runge_kutta third order*/
		static matrix<T> rk3(const matrix<T>&,matrix<T> (*)(const matrix<T>&),T,T,unsigned int);
		static matrix<T> rk3(const matrix<T>&,matrix<T> (*)(const matrix<T>&,T),T,T,unsigned int);
		static matrix<T> rk3(const matrix<T>&,matrix<T> (*)(const matrix<T>&),T,T,char [],unsigned int);
		static matrix<T> rk3(const matrix<T>&,matrix<T> (*)(const matrix<T>&,T),T,T,char [],unsigned int);
		
		/*Runge_kutta fourth order*/
		static matrix<T> rk4(const matrix<T>&,matrix<T> (*)(const matrix<T>&),T,T,unsigned int);
		static matrix<T> rk4(const matrix<T>&,matrix<T> (*)(const matrix<T>&,T),T,T,unsigned int);
		static matrix<T> rk4(const matrix<T>&,matrix<T> (*)(const matrix<T>&),T,T,char [],unsigned int);
		static matrix<T> rk4(const matrix<T>&,matrix<T> (*)(const matrix<T>&,T),T,T,char [],unsigned int);
		
		/*Implict Euler*/
		static matrix<T> ieuler(const matrix<T>&,matrix<T> (*)(const matrix<T>&),T,T,unsigned int,T,T,unsigned int);
		static matrix<T> ieuler(const matrix<T>&,matrix<T> (*)(const matrix<T>&,T),T,T,unsigned int,T,T,unsigned int);
		static matrix<T> ieuler(const matrix<T>&,matrix<T> (*)(const matrix<T>&),T,T,char [],unsigned int,T,T,unsigned int);
		static matrix<T> ieuler(const matrix<T>&,matrix<T> (*)(const matrix<T>&,T),T,T,char [],unsigned int,T,T,unsigned int);
		
		/*Crank Nicolson*/
		static matrix<T> crank_nicolson(const matrix<T>&,matrix<T> (*)(const matrix<T>&),T,T,unsigned int,T,T,unsigned int);
		static matrix<T> crank_nicolson(const matrix<T>&,matrix<T> (*)(const matrix<T>&,T),T,T,unsigned int,T,T,unsigned int);
		static matrix<T> crank_nicolson(const matrix<T>&,matrix<T> (*)(const matrix<T>&),T,T,char [],unsigned int,T,T,unsigned int);
		static matrix<T> crank_nicolson(const matrix<T>&,matrix<T> (*)(const matrix<T>&,T),T,T,char [],unsigned int,T,T,unsigned int);
		
		/*Crank Nicolson tridiag*/
		static matrix<T> crank_nicolson_tridiag(const matrix<T>&,matrix<T> (*)(const matrix<T>&),T,T,unsigned int,T,T,unsigned int);
		static matrix<T> crank_nicolson_tridiag(const matrix<T>&,matrix<T> (*)(const matrix<T>&,T),T,T,unsigned int,T,T,unsigned int);
		static matrix<T> crank_nicolson_tridiag(const matrix<T>&,matrix<T> (*)(const matrix<T>&),T,T,char [],unsigned int,T,T,unsigned int);
		static matrix<T> crank_nicolson_tridiag(const matrix<T>&,matrix<T> (*)(const matrix<T>&,T),T,T,char [],unsigned int,T,T,unsigned int);
		
		/*Appoximate for tridiag matrix ~ usually second dervitive pdes*/
		static matrix<T> ieuler_tridiag(const matrix<T>&,matrix<T> (*)(const matrix<T>&),T,T,unsigned int,T,T,unsigned int);
		static matrix<T> ieuler_tridiag(const matrix<T>&,matrix<T> (*)(const matrix<T>&,T),T,T,unsigned int,T,T,unsigned int);
		static matrix<T> ieuler_tridiag(const matrix<T>&,matrix<T> (*)(const matrix<T>&),T,T,char [],unsigned int,T,T,unsigned int);
		static matrix<T> ieuler_tridiag(const matrix<T>&,matrix<T> (*)(const matrix<T>&,T),T,T,char [],unsigned int,T,T,unsigned int);
		
		/*Implict Euler for when the Jacobian is given functions to come soon*/
	};

	/*Euler*/
	template<class T>
	matrix<T> ode<T>::euler(const matrix<T> &init,matrix<T> (*F)(const matrix<T>&,T),T end_time,T dt,unsigned int samples) //Techincally starts at t = 0
	{
		/*Checks to make sure the init is a vector*/
		if(init.cols() != 1) {std::cout << "Runtime Error\nFunction must be a vector\n"; exit(1);} else {}
		/*Checks to see if we can get the correct number of samples in how many steps*/
		if(samples > (unsigned int)((T)end_time/dt)) {std::cout << "Runtime Error\nNumber of samples is greater then number of runs\n";exit(1);} else {}
		/*Copy the inital condtion*/
		matrix<T> inital = init;
		
		/*Set up the output file*/
		FILE *file = fopen("output.txt","w");
		if(!file) {std::cout << "Cannot open file\n";exit(1);} else {}
		
		/*Counter for when to record the data*/
		unsigned int steps = (unsigned int) (end_time/dt)/(T)samples;
		
		/*Start solving*/
		unsigned int i = 0;
		while(dt*(i) <= end_time)
		{	
			if(steps == (unsigned int) (end_time/dt)/(T)samples)
			{
				/*Tells how far into the program it is*/
				std::cout << 100*(dt*(T)i/(T)end_time) << "% Done\n";
				
				/*Print results to text file*/
				fprintf(file,"%.24f\t",dt*i);
				for(int i = 0; i < inital.numel(); i++)
				{
					fprintf(file,"%.24f\t",inital(i));
				}
				fprintf(file,"\n");
				steps = 0;
			}
			else {steps++;}
			
			/*Solves for the next time step*/
			inital += dt*F(inital,dt*i);
			/*Move to the next time step*/
			i++;
		}
		/*Close the File and return the final solution*/
		fclose(file);
		std::cout << 100 << "% Done\n";
		return inital;
	}
	
	template<class T>
	matrix<T> ode<T>::euler(const matrix<T> &init,matrix<T> (*F)(const matrix<T>&),T end_time,T dt,unsigned int samples) //Techincally starts at t = 0
	{
		/*Checks to make sure the init is a vector*/
		if(init.cols() != 1) {std::cout << "Runtime Error\nFunction must be a vector\n"; exit(1);} else {}
		/*Checks to see if we can get the correct number of samples in how many steps*/
		if(samples > (unsigned int)((T)end_time/dt)) {std::cout << "Runtime Error\nNumber of samples is greater then number of runs\n";exit(1);} else {}
		/*Copy the inital condtion*/
		matrix<T> inital = init;
		
		/*Set up the output file*/
		FILE *file = fopen("output.txt","w");
		if(!file) {std::cout << "Cannot open file\n";exit(1);} else {}
		
		/*Counter for when to record the data*/
		unsigned int steps = (unsigned int) (end_time/dt)/(T)samples;
		
		/*Start solving*/
		unsigned int i = 0;
		while(dt*(i) <= end_time)
		{
			if(steps == (unsigned int) (end_time/dt)/(T)samples)
			{
				/*Tells how far into the program it is*/
				std::cout << 100*(dt*(double)i/(double)end_time) << "% Done\n";
				
				/*Print results to text file*/
				fprintf(file,"%.24f\t",dt*i);
				for(int i = 0; i < inital.numel(); i++)
				{
					fprintf(file,"%.24f\t",inital(i));
				}
				fprintf(file,"\n");
				steps = 0;
			}
			else {steps++;}
			
			/*Solves for the next time step*/
			inital += dt*F(inital);
			
			/*Move to the next time step*/
			i++;
			
		}
		/*Close the File and return the final solution*/
		fclose(file);
		std::cout << 100 << "% Done\n";
		return inital;
	}

	template<class T>
	matrix<T> ode<T>::euler(const matrix<T> &init,matrix<T> (*F)(const matrix<T>&,T),T end_time,T dt,char name[],unsigned int samples) //Techincally starts at t = 0
	{
		/*Checks to make sure the init is a vector*/
		if(init.cols() != 1) {std::cout << "Runtime Error\nFunction must be a vector\n"; exit(1);} else {}
		/*Checks to see if we can get the correct number of samples in how many steps*/
		if(samples > (unsigned int)((T)end_time/dt)) {std::cout << "Runtime Error\nNumber of samples is greater then number of runs\n";exit(1);} else {}
		/*Copy the inital condtion*/
		matrix<T> inital = init;
		
		/*Set up the output file*/
		FILE *file = fopen(name,"w");
		if(!file) {std::cout << "Cannot open file\n";exit(1);} else {}
		
		/*Counter for when to record the data*/
		unsigned int steps = (unsigned int) (end_time/dt)/(T)samples;
		
		/*Start solving*/
		unsigned int i = 0;
		while(dt*(i) <= end_time)
		{
			if(steps == (unsigned int) (end_time/dt)/(T)samples)
			{
				/*Tells how far into the program it is*/
				std::cout << 100*(dt*(double)i/(double)end_time) << "% Done\n";
				
				/*Print results to text file*/
				fprintf(file,"%.24f\t",dt*i);
				for(int i = 0; i < inital.numel(); i++)
				{
					fprintf(file,"%.24f\t",inital(i));
				}
				fprintf(file,"\n");
				steps = 0;
			}
			else {steps++;}
			
			/*Solves for the next time step*/
			inital += dt*F(inital,dt*i);
			
			/*Move to the next time step*/
			i++;
		}
		/*Close the File and return the final solution*/
		fclose(file);
		std::cout << 100 << "% Done\n";
		return inital;
	}
	
	template<class T>
	matrix<T> ode<T>::euler(const matrix<T> &init,matrix<T> (*F)(const matrix<T>&),T end_time,T dt,char name[],unsigned int samples) //Techincally starts at t = 0
	{
		/*Checks to make sure the init is a vector*/
		if(init.cols() != 1) {std::cout << "Runtime Error\nFunction must be a vector\n"; exit(1);} else {}
		/*Checks to see if we can get the correct number of samples in how many steps*/
		if(samples > (unsigned int)((T)end_time/dt)) {std::cout << "Runtime Error\nNumber of samples is greater then number of runs\n";exit(1);} else {}
		/*Copy the inital condtion*/
		matrix<T> inital = init;
		
		/*Set up the output file*/
		FILE *file = fopen(name,"w");
		if(!file) {std::cout << "Cannot open file\n";exit(1);} else {}
		
		/*Counter for when to record the data*/
		unsigned int steps = (unsigned int) (end_time/dt)/(T)samples;
		
		/*Start solving*/
		unsigned int i = 0;
		while(dt*(i) <= end_time)
		{
			if(steps == (unsigned int) (end_time/dt)/(T)samples)
			{
				/*Tells how far into the program it is*/
				std::cout << 100*(dt*(double)i/(double)end_time) << "% Done\n";
				
				/*Print results to text file*/
				fprintf(file,"%.24f\t",dt*i);
				for(int i = 0; i < inital.numel(); i++)
				{
					fprintf(file,"%.24f\t",inital(i));
				}
				fprintf(file,"\n");
				steps = 0;
			}
			else {steps++;}
			
			/*Solves for the next time step*/
			inital += dt*F(inital);
			
			/*Move to the next time step*/
			i++;
			
		}
		/*Close the File and return the final solution*/
		fclose(file);
		std::cout << 100 << "% Done\n";
		return inital;
	}

	/*runge kutta 2*/
	
	template<class T>
	matrix<T> ode<T>::rk2(const matrix<T> &init,matrix<T> (*F)(const matrix<T>&,T),T end_time,T dt,unsigned int samples)
	{
		/*Checks to make sure the init is a vector*/
		if(init.cols() != 1) {std::cout << "Runtime Error\nFunction must be a vector\n"; exit(1);} else {}
		/*Checks to see if we can get the correct number of samples in how many steps*/
		if(samples > (unsigned int)((T)end_time/dt)) {std::cout << "Runtime Error\nNumber of samples is greater then number of runs\n";exit(1);} else {}
		/*Copy the inital condtion*/
		matrix<T> inital = init;
		
		/*Set up the output file*/
		FILE *file = fopen("output.txt","w");
		if(!file) {std::cout << "Cannot open file\n";exit(1);} else {}
		
		/*Counter for when to record the data*/
		unsigned int steps = (unsigned int) (end_time/dt)/(T)samples;
		
		/*Start solving*/
		unsigned int i = 0;
		while(dt*(i) <= end_time)
		{
			if(steps == (unsigned int) (end_time/dt)/(T)samples)
			{
				/*Tells how far into the program it is*/
				std::cout << 100*(dt*(double)i/(double)end_time) << "% Done\n";
				
				/*Print results to text file*/
				fprintf(file,"%.24f\t",dt*i);
				for(int i = 0; i < inital.numel(); i++)
				{
					fprintf(file,"%.24f\t",inital(i));
				}
				fprintf(file,"\n");
				steps = 0;
			}
			else {steps++;}
			
			/*Solves for the next time step*/
			matrix<T> k1 = inital + .5*dt*F(inital,dt*i);
			inital += dt*F(k1,dt*i + .5*dt);
			
			/*Update time steps*/
			i++;
		}
		/*Close the File and return the final solution*/
		fclose(file);
		std::cout << 100 << "% Done\n";
		return inital;
	}

	template<class T>
	matrix<T> ode<T>::rk2(const matrix<T> &init,matrix<T> (*F)(const matrix<T>&),T end_time,T dt,unsigned int samples)
	{
		/*Checks to make sure the init is a vector*/
		if(init.cols() != 1) {std::cout << "Runtime Error\nFunction must be a vector\n"; exit(1);} else {}
		/*Checks to see if we can get the correct number of samples in how many steps*/
		if(samples > (unsigned int)((T)end_time/dt)) {std::cout << "Runtime Error\nNumber of samples is greater then number of runs\n";exit(1);} else {}
		/*Copy the inital condtion*/
		matrix<T> inital = init;
		
		/*Set up the output file*/
		FILE *file = fopen("output.txt","w");
		if(!file) {std::cout << "Cannot open file\n";exit(1);} else {}
		
		/*Counter for when to record the data*/
		unsigned int steps = (unsigned int) (end_time/dt)/(T)samples;
		
		/*Start solving*/
		unsigned int i = 0;
		while(dt*(i) <= end_time)
		{
			if(steps == (unsigned int) (end_time/dt)/(T)samples)
			{
				/*Tells how far into the program it is*/
				std::cout << 100*(dt*(double)i/(double)end_time) << "% Done\n";
				
				/*Print results to text file*/
				fprintf(file,"%.24f\t",dt*i);
				for(int i = 0; i < inital.numel(); i++)
				{
					fprintf(file,"%.24f\t",inital(i));
				}
				fprintf(file,"\n");
				steps = 0;
			}
			else {steps++;}
			
			/*Solves for the next time step*/
			matrix<T> k1 = inital + .5*dt*F(inital);
			inital += dt*F(k1);
			
			/*Update time steps*/
			i++;
		}
		/*Close the File and return the final solution*/
		fclose(file);
		std::cout << 100 << "% Done\n";
		return inital;
	}

	template<class T>
	matrix<T> ode<T>::rk2(const matrix<T> &init,matrix<T> (*F)(const matrix<T>&,T),T end_time,T dt,char name[],unsigned int samples)
	{
		/*Checks to make sure the init is a vector*/
		if(init.cols() != 1) {std::cout << "Runtime Error\nFunction must be a vector\n"; exit(1);} else {}
		/*Checks to see if we can get the correct number of samples in how many steps*/
		if(samples > (unsigned int)((T)end_time/dt)) {std::cout << "Runtime Error\nNumber of samples is greater then number of runs\n";exit(1);} else {}
		/*Copy the inital condtion*/
		matrix<T> inital = init;
		
		/*Set up the output file*/
		FILE *file = fopen(name,"w");
		if(!file) {std::cout << "Cannot open file\n";exit(1);} else {}
		
		/*Counter for when to record the data*/
		unsigned int steps = (unsigned int) (end_time/dt)/(T)samples;
		
		/*Start solving*/
		unsigned int i = 0;
		while(dt*(i) <= end_time)
		{
			if(steps == (unsigned int) (end_time/dt)/(T)samples)
			{
				/*Tells how far into the program it is*/
				std::cout << 100*(dt*(double)i/(double)end_time) << "% Done\n";
				
				/*Print results to text file*/
				fprintf(file,"%.24f\t",dt*i);
				for(int i = 0; i < inital.numel(); i++)
				{
					fprintf(file,"%.24f\t",inital(i));
				}
				fprintf(file,"\n");
				steps = 0;
			}
			else {steps++;}
			
			/*Solves for the next time step*/
			matrix<T> k1 = inital + .5*dt*F(inital,dt*i);
			inital += dt*F(k1,dt*i + .5*dt);
			
			/*Update time steps*/
			i++;
		}
		/*Close the File and return the final solution*/
		fclose(file);
		std::cout << 100 << "% Done\n";
		return inital;
	}

	template<class T>
	matrix<T> ode<T>::rk2(const matrix<T> &init,matrix<T> (*F)(const matrix<T>&),T end_time,T dt,char name[],unsigned int samples)
	{
		/*Checks to make sure the init is a vector*/
		if(init.cols() != 1) {std::cout << "Runtime Error\nFunction must be a vector\n"; exit(1);} else {}
		/*Checks to see if we can get the correct number of samples in how many steps*/
		if(samples > (unsigned int)((T)end_time/dt)) {std::cout << "Runtime Error\nNumber of samples is greater then number of runs\n";exit(1);} else {}
		/*Copy the inital condtion*/
		matrix<T> inital = init;
		
		/*Set up the output file*/
		FILE *file = fopen(name,"w");
		if(!file) {std::cout << "Cannot open file\n";exit(1);} else {}
		
		/*Counter for when to record the data*/
		unsigned int steps = (unsigned int) (end_time/dt)/(T)samples;
		
		/*Start solving*/
		unsigned int i = 0;
		while(dt*(i) <= end_time)
		{
			if(steps == (unsigned int) (end_time/dt)/(T)samples)
			{
				/*Tells how far into the program it is*/
				std::cout << 100*(dt*(double)i/(double)end_time) << "% Done\n";
				
				/*Print results to text file*/
				fprintf(file,"%.24f\t",dt*i);
				for(int i = 0; i < inital.numel(); i++)
				{
					fprintf(file,"%.24f\t",inital(i));
				}
				fprintf(file,"\n");
				steps = 0;
			}
			else {steps++;}
			
			/*Solves for the next time step*/
			matrix<T> k1 = inital + .5*dt*F(inital);
			inital += dt*F(k1);

			/*Update time steps*/
			i++;
		}
		/*Close the File and return the final solution*/
		fclose(file);
		std::cout << 100 << "% Done\n";
		return inital;
	}

	/*Runge Kutta third order*/
	
	template<class T>
	matrix<T> ode<T>::rk3(const matrix<T> &init,matrix<T> (*F)(const matrix<T>&),T end_time,T dt,unsigned int samples)
	{
		/*Checks to make sure the init is a vector*/
		if(init.cols() != 1) {std::cout << "Runtime Error\nFunction must be a vector\n"; exit(1);} else {}
		/*Checks to see if we can get the correct number of samples in how many steps*/
		if(samples > (unsigned int)((T)end_time/dt)) {std::cout << "Runtime Error\nNumber of samples is greater then number of runs\n";exit(1);} else {}
		/*Copy the inital condtion*/
		matrix<T> inital = init;
		
		/*Set up the output file*/
		FILE *file = fopen("output.txt","w");
		if(!file) {std::cout << "Cannot open file\n";exit(1);} else {}
		
		/*Counter for when to record the data*/
		unsigned int steps = (unsigned int) (end_time/dt)/(T)samples;
		
		/*Start solving*/
		unsigned int i = 0;
		while(dt*(i) <= end_time)
		{
			if(steps == (unsigned int) (end_time/dt)/(T)samples)
			{
				/*Tells how far into the program it is*/
				std::cout << 100*(dt*(double)i/(double)end_time) << "% Done\n";
				
				/*Print results to text file*/
				fprintf(file,"%.24f\t",dt*i);
				for(int i = 0; i < inital.numel(); i++)
				{
					fprintf(file,"%.24f\t",inital(i));
				}
				fprintf(file,"\n");
				steps = 0;
			}
			else {steps++;}
			
			/*Solves for the next time step*/
			matrix<T> k1 = dt*F(inital);
			matrix<T> k2 = dt*F(inital + k1/(T)2);
			matrix<T> k3 = dt*F(inital - k1 + (T)2*k2);
			inital += (T)(1./6.)*(k1 + (T)4*k2 + k3);
			
			/*Update time step*/
			i++;
		}
		/*Close the File and return the final solution*/
		fclose(file);
		std::cout << 100 << "% Done\n";
		return inital;
	}

	template<class T>
	matrix<T> ode<T>::rk3(const matrix<T> &init,matrix<T> (*F)(const matrix<T>&,T),T end_time,T dt,unsigned int samples)
	{
		/*Checks to make sure the init is a vector*/
		if(init.cols() != 1) {std::cout << "Runtime Error\nFunction must be a vector\n"; exit(1);} else {}
		/*Checks to see if we can get the correct number of samples in how many steps*/
		if(samples > (unsigned int)((T)end_time/dt)) {std::cout << "Runtime Error\nNumber of samples is greater then number of runs\n";exit(1);} else {}
		/*Copy the inital condtion*/
		matrix<T> inital = init;
		
		/*Set up the output file*/
		FILE *file = fopen("output.txt","w");
		if(!file) {std::cout << "Cannot open file\n";exit(1);} else {}
		
		/*Counter for when to record the data*/
		unsigned int steps = (unsigned int) (end_time/dt)/(T)samples;
		
		/*Start solving*/
		unsigned int i = 0;
		while(dt*(i) <= end_time)
		{
			if(steps == (unsigned int) (end_time/dt)/(T)samples)
			{
				/*Tells how far into the program it is*/
				std::cout << 100*(dt*(double)i/(double)end_time) << "% Done\n";
				
				/*Print results to text file*/
				fprintf(file,"%.24f\t",dt*i);
				for(int i = 0; i < inital.numel(); i++)
				{
					fprintf(file,"%.24f\t",inital(i));
				}
				fprintf(file,"\n");
				steps = 0;
			}
			else {steps++;}
			
			/*Solves for the next time step*/
			matrix<T> k1 = dt*F(inital,dt*i);
			matrix<T> k2 = dt*F(inital + k1/(T)2,dt*i + dt/(T)2);
			matrix<T> k3 = dt*F(inital - k1 + (T)2*k2,dt*i + dt);
			inital += (T)(1./6.)*(k1 + (T)4*k2 + k3);

			/*Update time step*/
			i++;
		}
		/*Close the File and return the final solution*/
		fclose(file);
		std::cout << 100 << "% Done\n";
		return inital;
	}

	template<class T>
	matrix<T> ode<T>::rk3(const matrix<T> &init,matrix<T> (*F)(const matrix<T>&),T end_time,T dt,char name[],unsigned int samples)
	{
		/*Checks to make sure the init is a vector*/
		if(init.cols() != 1) {std::cout << "Runtime Error\nFunction must be a vector\n"; exit(1);} else {}
		/*Checks to see if we can get the correct number of samples in how many steps*/
		if(samples > (unsigned int)((T)end_time/dt)) {std::cout << "Runtime Error\nNumber of samples is greater then number of runs\n";exit(1);} else {}
		/*Copy the inital condtion*/
		matrix<T> inital = init;
		
		/*Set up the output file*/
		FILE *file = fopen(name,"w");
		if(!file) {std::cout << "Cannot open file\n";exit(1);} else {}
		
		/*Counter for when to record the data*/
		unsigned int steps = (unsigned int) (end_time/dt)/(T)samples;
		
		/*Start solving*/
		unsigned int i = 0;
		while(dt*(i) <= end_time)
		{
			if(steps ==(unsigned int) (end_time/dt)/(T)samples)
			{
				/*Tells how far into the program it is*/
				std::cout << 100*(dt*(double)i/(double)end_time) << "% Done\n";
				
				/*Print results to text file*/
				fprintf(file,"%.24f\t",dt*i);
				for(int i = 0; i < inital.numel(); i++)
				{
					fprintf(file,"%.24f\t",inital(i));
				}
				fprintf(file,"\n");
				steps = 0;
			}
			else {steps++;}
			
			/*Solves for the next time step*/
			matrix<T> k1 = dt*F(inital);
			matrix<T> k2 = dt*F(inital + k1/(T)2);
			matrix<T> k3 = dt*F(inital - k1 + (T)2*k2);
			inital += (T)(1./6.)*(k1 + (T)4*k2 + k3);

			/*Update time step*/
			i++;
		}
		/*Close the File and return the final solution*/
		fclose(file);
		std::cout << 100 << "% Done\n";
		return inital;
	}

	template<class T>
	matrix<T> ode<T>::rk3(const matrix<T> &init,matrix<T> (*F)(const matrix<T>&,T),T end_time,T dt,char name[],unsigned int samples)
	{
		/*Checks to make sure the init is a vector*/
		if(init.cols() != 1) {std::cout << "Runtime Error\nFunction must be a vector\n"; exit(1);} else {}
		/*Checks to see if we can get the correct number of samples in how many steps*/
		if(samples > (unsigned int)((T)end_time/dt)) {std::cout << "Runtime Error\nNumber of samples is greater then number of runs\n";exit(1);} else {}
		/*Copy the inital condtion*/
		matrix<T> inital = init;
		
		/*Set up the output file*/
		FILE *file = fopen(name,"w");
		if(!file) {std::cout << "Cannot open file\n";exit(1);} else {}
		
		/*Counter for when to record the data*/
		unsigned int steps = (unsigned int) (end_time/dt)/(T)samples;
		
		/*Start solving*/
		unsigned int i = 0;
		while(dt*(i) <= end_time)
		{
			if(steps == (unsigned int) (end_time/dt)/(T)samples)
			{
				/*Tells how far into the program it is*/
				std::cout << 100*(dt*(double)i/(double)end_time) << "% Done\n";
				
				/*Print results to text file*/
				fprintf(file,"%.24f\t",dt*i);
				for(int i = 0; i < inital.numel(); i++)
				{
					fprintf(file,"%.24f\t",inital(i));
				}
				fprintf(file,"\n");
				steps = 0;
			}
			else {steps++;}
			
			/*Solves for the next time step*/
			matrix<T> k1 = dt*F(inital,dt*i);
			matrix<T> k2 = dt*F(inital + k1/(T)2,dt*i + dt/(T)2);
			matrix<T> k3 = dt*F(inital - k1 + (T)2*k2,dt*i + dt);
			inital += (T)(1./6.)*(k1 + (T)4*k2 + k3);
			
			/*Update time step*/
			i++;
		}
		/*Close the File and return the final solution*/
		fclose(file);
		std::cout << 100 << "% Done\n";
		return inital;
	}

	/*Runge kutta fourth order methods*/
	
	template<class T>
	matrix<T> ode<T>::rk4(const matrix<T> &init,matrix<T> (*F)(const matrix<T>&),T end_time,T dt,unsigned int samples)
	{
		/*Checks to make sure the init is a vector*/
		if(init.cols() != 1) {std::cout << "Runtime Error\nFunction must be a vector\n"; exit(1);} else {}
		/*Checks to see if we can get the correct number of samples in how many steps*/
		if(samples > (unsigned int)((T)end_time/dt)) {std::cout << "Runtime Error\nNumber of samples is greater then number of runs\n";exit(1);} else {}
		/*Copy the inital condtion*/
		matrix<T> inital = init;
		
		/*Initalize the k_is*/
		matrix<T> k1(init.rows(),1);
		matrix<T> k2(init.rows(),1);
		matrix<T> k3(init.rows(),1);
		matrix<T> k4(init.rows(),1);
		
		/*Set up the output file*/
		FILE *file = fopen("output.txt","w");
		if(!file) {std::cout << "Cannot open file\n";exit(1);} else {}
		
		/*Counter for when to record the data*/
		unsigned int steps = (unsigned int) (end_time/dt)/(T)samples;
		
		/*Start solving*/
		unsigned int i = 0;
		while(dt*(i) <= end_time)
		{
			if(steps == (unsigned int) (end_time/dt)/(T)samples)
			{
				/*Tells how far into the program it is*/
				std::cout << 100*(dt*(double)i/(double)end_time) << "% Done\n";
				
				/*Print results to text file*/
				fprintf(file,"%.24f\t",dt*i);
				for(int i = 0; i < inital.numel(); i++)
				{
					fprintf(file,"%.24f\t",inital(i));
				}
				fprintf(file,"\n");
				steps = 0;
			}
			else {steps++;}
			
			/*Solves for the next time step*/
			k1 = F(inital);
			k2 = F(inital + (dt/(T)2)*k1);
			k3 = F(inital + (dt/(T)2)*k2);
			k4 = F(inital + dt*k3);
			inital += dt*(T)(1./6.)*(k1 + (T)2*k2 + (T)2*k3 + k4);

			/*Update to next time step*/
			i++;
		}
		/*Close the File and return the final solution*/
		fclose(file);
		std::cout << 100 << "% Done\n";
		return inital;
	}

	template<class T>
	matrix<T> ode<T>::rk4(const matrix<T> &init,matrix<T> (*F)(const matrix<T>&,T),T end_time,T dt,unsigned int samples)
	{
		/*Checks to make sure the init is a vector*/
		if(init.cols() != 1) {std::cout << "Runtime Error\nFunction must be a vector\n"; exit(1);} else {}
		/*Checks to see if we can get the correct number of samples in how many steps*/
		if(samples > (unsigned int)((T)end_time/dt)) {std::cout << "Runtime Error\nNumber of samples is greater then number of runs\n";exit(1);} else {}
		/*Copy the inital condtion*/
		matrix<T> inital = init;
		
		/*Initalize the k_is*/
		matrix<T> k1(init.rows(),1);
		matrix<T> k2(init.rows(),1);
		matrix<T> k3(init.rows(),1);
		matrix<T> k4(init.rows(),1);
		
		/*Set up the output file*/
		FILE *file = fopen("output.txt","w");
		if(!file) {std::cout << "Cannot open file\n";exit(1);} else {}
		
		/*Counter for when to record the data*/
		unsigned int steps = (unsigned int) (end_time/dt)/(T)samples;
		
		/*Start solving*/
		unsigned int i = 0;
		while(dt*(i) <= end_time)
		{
			if(steps == (unsigned int) (end_time/dt)/(T)samples)
			{
				/*Tells how far into the program it is*/
				std::cout << 100*(dt*(double)i/(double)end_time) << "% Done\n";
				
				/*Print results to text file*/
				fprintf(file,"%.24f\t",dt*i);
				for(int i = 0; i < inital.numel(); i++)
				{
					fprintf(file,"%.24f\t",inital(i));
				}
				fprintf(file,"\n");
				steps = 0;
			}
			else {steps++;}
			
			/*Solves for the next time step*/
			k1 = F(inital,dt*i);
			k2 = F(inital + (dt/(T)2)*k1,dt*i + dt/(T)2);
			k3 = F(inital + (dt/(T)2)*k2,dt*i + dt/(T)2);
			k4 = F(inital + dt*k3,dt*i + dt);
			inital += dt*(T)(1./6.)*(k1 + (T)2*k2 + (T)2*k3 + k4);

			/*Update to next time step*/
			i++;
		}
		/*Close the File and return the final solution*/
		fclose(file);
		std::cout << 100 << "% Done\n";
		return inital;
	}

	template<class T>
	matrix<T> ode<T>::rk4(const matrix<T> &init,matrix<T> (*F)(const matrix<T>&),T end_time,T dt,char name[],unsigned int samples)
	{
		/*Checks to make sure the init is a vector*/
		if(init.cols() != 1) {std::cout << "Runtime Error\nFunction must be a vector\n"; exit(1);} else {}
		/*Checks to see if we can get the correct number of samples in how many steps*/
		if(samples > (unsigned int)((T)end_time/dt)) {std::cout << "Runtime Error\nNumber of samples is greater then number of runs\n";exit(1);} else {}
		/*Copy the inital condtion*/
		matrix<T> inital = init;
		
		/*Initalize the k_is*/
		matrix<T> k1(init.rows(),1);
		matrix<T> k2(init.rows(),1);
		matrix<T> k3(init.rows(),1);
		matrix<T> k4(init.rows(),1);
		
		/*Set up the output file*/
		FILE *file = fopen(name,"w");
		if(!file) {std::cout << "Cannot open file\n";exit(1);} else {}
		
		/*Counter for when to record the data*/
		unsigned int steps = (unsigned int) (end_time/dt)/(T)samples;
		
		/*Start solving*/
		unsigned int i = 0;
		while(dt*(i) <= end_time)
		{
			if(steps == (unsigned int) (end_time/dt)/(T)samples)
			{
				/*Tells how far into the program it is*/
				std::cout << 100*(dt*(double)i/(double)end_time) << "% Done\n";
				
				/*Print results to text file*/
				fprintf(file,"%.24f\t",dt*i);
				for(int i = 0; i < inital.numel(); i++)
				{
					fprintf(file,"%.24f\t",inital(i));
				}
				fprintf(file,"\n");
				steps = 0;
			}
			else {steps++;}
			
			/*Solves for the next time step*/
			k1 = F(inital);
			k2 = F(inital + (dt/(T)2)*k1);
			k3 = F(inital + (dt/(T)2)*k2);
			k4 = F(inital + dt*k3);
			inital += dt*(T)(1./6.)*(k1 + (T)2*k2 + (T)2*k3 + k4);

			/*Update to next time step*/
			i++;
		}
		/*Close the File and return the final solution*/
		fclose(file);
		std::cout << 100 << "% Done\n";
		return inital;
	}

	template<class T>
	matrix<T> ode<T>::rk4(const matrix<T> &init,matrix<T> (*F)(const matrix<T>&,T),T end_time,T dt,char name[],unsigned int samples)
	{
		/*Checks to make sure the init is a vector*/
		if(init.cols() != 1) {std::cout << "Runtime Error\nFunction must be a vector\n"; exit(1);} else {}
		/*Checks to see if we can get the correct number of samples in how many steps*/
		if(samples > (unsigned int)((T)end_time/dt)) {std::cout << "Runtime Error\nNumber of samples is greater then number of runs\n";exit(1);} else {}
		/*Copy the inital condtion*/
		matrix<T> inital = init;
		
		/*Initalize the k_is*/
		matrix<T> k1(init.rows(),1);
		matrix<T> k2(init.rows(),1);
		matrix<T> k3(init.rows(),1);
		matrix<T> k4(init.rows(),1);

		/*Set up the output file*/
		FILE *file = fopen(name,"w");
		if(!file) {std::cout << "Cannot open file\n";exit(1);} else {}
		
		/*Counter for when to record the data*/
		unsigned int steps = (unsigned int) (end_time/dt)/(T)samples;
		
		/*Start solving*/
		unsigned int i = 0;
		while(dt*(i) <= end_time)
		{
			if(steps == (unsigned int) (end_time/dt)/(T)samples)
			{
				/*Tells how far into the program it is*/
				std::cout << 100*(dt*(double)i/(double)end_time) << "% Done\n";
				
				/*Print results to text file*/
				fprintf(file,"%.24f\t",dt*i);
				for(int i = 0; i < inital.numel(); i++)
				{
					fprintf(file,"%.24f\t",inital(i));
				}
				fprintf(file,"\n");
				steps = 0;
			}
			else {steps++;}
			
			/*Solves for the next time step*/
			k1 = F(inital,dt*i);
			k2 = F(inital + (dt/(T)2)*k1,dt*i + dt/(T)2);
			k3 = F(inital + (dt/(T)2)*k2,dt*i + dt/(T)2);
			k4 = F(inital + dt*k3,dt*i + dt);
			inital += dt*(T)(1./6.)*(k1 + (T)2*k2 + (T)2*k3 + k4);

			/*Update to next time step*/
			i++;
		}
		/*Close the File and return the final solution*/
		fclose(file);
		std::cout << 100 << "% Done\n";
		return inital;
	}

	/*Implict Methods*/
	/*Must set the form in proper form where it equals a residue*/
	template<class T>
	matrix<T> ode<T>::ieuler(const matrix<T> &init,matrix<T> (*F)(const matrix<T>&),T end_time,T dt,unsigned int samples,T tol,T step_size,unsigned int max_iter)
	{
		/*Checks to make sure the init is a vector*/
		if(init.cols() != 1) {std::cout << "Runtime Error\nFunction must be a vector\n"; exit(1);} else {}
		/*Checks to see if we can get the correct number of samples in how many steps*/
		if(samples > (unsigned int)((T)end_time/dt)) {std::cout << "Runtime Error\nNumber of samples is greater then number of runs\n";exit(1);} else {}
		/*Copy the inital condtion*/
		matrix<T> inital = init;
		matrix<T> J = matrix<T>::zeros(init.numel(),init.numel());
		matrix<T> fv = matrix<T>::zeros(init.numel(),1);
	
		/*Get the size of the system*/
		const int N = inital.numel();
	
		/*Set up the output file*/
		FILE *file = fopen("output.txt","w");
		if(!file) {std::cout << "Cannot open file\n"; exit(1);} else {}
	
		/*Counter for when to record the data*/
		unsigned int steps = (unsigned int) (end_time/dt)/(T)samples;
	
		/*Start solving*/
		unsigned int i = 0;
		while(dt*(i) <= end_time)
		{
			if(steps == (unsigned int) (end_time/dt)/(T)samples)
			{
				/*Tells how far into the program it is*/
				std::cout << 100*(dt*(double)i/(double)end_time) << "% Done\n";
			
				/*Print results to text file*/
				fprintf(file,"%.14f\t",dt*i);
				for(int i = 0; i < inital.numel(); i++)
				{
					fprintf(file,"%.14f\t",inital(i));
				}
				fprintf(file,"\n");
				steps = 0;
			}
			else {steps++;}
		
			/*Solves for the next time step*/
			/*Use newtons root finding method to solve for the future timestep*/
		
			/*Our guess values*/
			matrix<T> g1 = inital + 100.*matrix<T>::ones(inital.numel(),1);
			matrix<T> g2 = inital;
			/*Start the root finding method*/
			unsigned int count = 0;
			while(norm_2(g1 - g2)/norm_2(g1) > tol && count++ < max_iter)
			{
				/*Update values*/
				g1 = g2;
				/*The Jacobian and function vector of the resiude*/
			
				/*Initalize the Jacobian*/
				//#pragma omp parallel for
				for(int j = 0; j < N; j++)
				{
					matrix<T> copy_inital = g1;
					copy_inital(j) += step_size;
					
					J.set_col(j,-dt*(F(copy_inital) - F(g1))/step_size);
				}
				J += matrix<T>::eye(N);
			
				/*Initalize the function vector*/
				fv = g1 - inital - dt*F(g1);
				g2 = g1 - gauss<T>::solve(J,fv);
				//std::cout << count << "\n";
			}
			if(norm_2(g2 - g1)/norm_2(g2) > tol) {std::cout << "Could not find solution in desired tolerance\n";} else{}
			inital = g2;
			
			/*Move to the next time step*/
			i++;
		}
		/*Close the File and return the final solution*/
		fclose(file);
		std::cout << 100 << "% Done\n";
		return inital;
	}
	
	template<class T>
	matrix<T> ode<T>::ieuler(const matrix<T> &init,matrix<T> (*F)(const matrix<T>&,T),T end_time,T dt,unsigned int samples,T tol,T step_size,unsigned int max_iter)
	{
		/*Checks to make sure the init is a vector*/
		if(init.cols() != 1) {std::cout << "Runtime Error\nFunction must be a vector\n"; exit(1);} else {}
		/*Checks to see if we can get the correct number of samples in how many steps*/
		if(samples > (unsigned int)((T)end_time/dt)) {std::cout << "Runtime Error\nNumber of samples is greater then number of runs\n";exit(1);} else {}
		/*Copy the inital condtion*/
		matrix<T> inital = init;
		matrix<T> J = matrix<T>::zeros(init.numel(),init.numel());
		matrix<T> fv = matrix<T>::zeros(init.numel(),1);
	
		/*Get the size of the system*/
		const int N = inital.numel();
	
		/*Set up the output file*/
		FILE *file = fopen("output.txt","w");
		if(!file) {std::cout << "Cannot open file\n"; exit(1);} else {}
	
		/*Counter for when to record the data*/
		unsigned int steps = (unsigned int) (end_time/dt)/(T)samples;
	
		/*Start solving*/
		unsigned int i = 0;
		while(dt*(i) <= end_time)
		{
			if(steps == (unsigned int) (end_time/dt)/(T)samples)
			{
				/*Tells how far into the program it is*/
				std::cout << 100*(dt*(double)i/(double)end_time) << "% Done\n";
			
				/*Print results to text file*/
				fprintf(file,"%.14f\t",dt*i);
				for(int i = 0; i < inital.numel(); i++)
				{
					fprintf(file,"%.14f\t",inital(i));
				}
				fprintf(file,"\n");
				steps = 0;
			}
			else {steps++;}
		
			/*Solves for the next time step*/
			/*Use newtons root finding method to solve for the future timestep*/
		
			/*Our guess values*/
			matrix<T> g1 = inital + 100.*matrix<T>::ones(inital.numel(),1);
			matrix<T> g2 = inital;
			/*Start the root finding method*/
			unsigned int count = 0;
			while(norm_2(g1 - g2)/norm_2(g1) > tol && count++ < max_iter)
			{
				/*Update values*/
				g1 = g2;
				/*The Jacobian and function vector of the resiude*/
			
				/*Initalize the Jacobian*/
				//#pragma omp parallel for
				for(int j = 0; j < N; j++)
				{
					matrix<T> copy_inital = g1;
					copy_inital(j) += step_size;
					
					J.set_col(j,-dt*(F(copy_inital,(T)dt*(i+1)) - F(g1,(T)dt*(i+1)))/step_size);
				}
				J += matrix<T>::eye(N);
			
				/*Initalize the function vector*/
				fv = g1 - inital - dt*F(g1,dt*(i+1));
				g2 = g1 - gauss<T>::solve(J,fv);
				//std::cout << count << "\n";
			}
			if(norm_2(g2 - g1)/norm_2(g2) > tol) {std::cout << "Could not find solution in desired tolerance\n";} else{}
			inital = g2;
			
			/*Move to the next time step*/
			i++;
		}
		/*Close the File and return the final solution*/
		fclose(file);
		std::cout << 100 << "% Done\n";
		return inital;
	}
	
	template<class T>
	matrix<T> ode<T>::ieuler(const matrix<T> &init,matrix<T> (*F)(const matrix<T>&),T end_time,T dt,char name[],unsigned int samples,T tol,T step_size,unsigned int max_iter)
	{
		/*Checks to make sure the init is a vector*/
		if(init.cols() != 1) {std::cout << "Runtime Error\nFunction must be a vector\n"; exit(1);} else {}
		/*Checks to see if we can get the correct number of samples in how many steps*/
		if(samples > (unsigned int)((T)end_time/dt)) {std::cout << "Runtime Error\nNumber of samples is greater then number of runs\n";exit(1);} else {}
		/*Copy the inital condtion*/
		matrix<T> inital = init;
		matrix<T> J = matrix<T>::zeros(init.numel(),init.numel());
		matrix<T> fv = matrix<T>::zeros(init.numel(),1);
		
		/*Get the size of the system*/
		const int N = inital.numel();
	
		/*Set up the output file*/
		FILE *file = fopen(name,"w");
		if(!file) {std::cout << "Cannot open file\n"; exit(1);} else {}
	
		/*Counter for when to record the data*/
		unsigned int steps = (unsigned int) (end_time/dt)/(T)samples;
	
		/*Start solving*/
		unsigned int i = 0;
		while(dt*(i) <= end_time)
		{
			if(steps == (unsigned int) (end_time/dt)/(T)samples)
			{
				/*Tells how far into the program it is*/
				std::cout << 100*(dt*(double)i/(double)end_time) << "% Done\n";
			
				/*Print results to text file*/
				fprintf(file,"%.14f\t",dt*i);
				for(int i = 0; i < inital.numel(); i++)
				{
					fprintf(file,"%.14f\t",inital(i));
				}
				fprintf(file,"\n");
				steps = 0;
			}
			else {steps++;}
		
			/*Solves for the next time step*/
			/*Use newtons root finding method to solve for the future timestep*/
		
			/*Our guess values*/
			matrix<T> g1 = inital + 100.*matrix<T>::ones(inital.numel(),1);
			matrix<T> g2 = inital;
			/*Start the root finding method*/
			unsigned int count = 0;
			while(norm_2(g1 - g2)/norm_2(g1) > tol && count++ < max_iter)
			{
				/*Update values*/
				g1 = g2;
				/*The Jacobian and function vector of the resiude*/
			
				/*Initalize the Jacobian*/
				//#pragma omp parallel for
				for(int j = 0; j < N; j++)
				{
					matrix<T> copy_inital = g1;
					copy_inital(j) += step_size;
					
					J.set_col(j,-dt*(F(copy_inital) - F(g1))/step_size);
				}
				J += matrix<T>::eye(N);
			
				/*Initalize the function vector*/
				fv = g1 - inital - dt*F(g1);
				g2 = g1 - gauss<T>::solve(J,fv);
				//std::cout << count << "\n";
			}
			if(norm_2(g2 - g1)/norm_2(g2) > tol) {std::cout << "Could not find solution in desired tolerance\n";} else{}
			inital = g2;
			
			/*Move to the next time step*/
			i++;
		}
		/*Close the File and return the final solution*/
		fclose(file);
		std::cout << 100 << "% Done\n";
		return inital;
	}
	
	template<class T>
	matrix<T> ode<T>::ieuler(const matrix<T> &init,matrix<T> (*F)(const matrix<T>&,T),T end_time,T dt,char name[],unsigned int samples,T tol,T step_size,unsigned int max_iter)
	{
		/*Checks to make sure the init is a vector*/
		if(init.cols() != 1) {std::cout << "Runtime Error\nFunction must be a vector\n"; exit(1);} else {}
		/*Checks to see if we can get the correct number of samples in how many steps*/
		if(samples > (unsigned int)((T)end_time/dt)) {std::cout << "Runtime Error\nNumber of samples is greater then number of runs\n";exit(1);} else {}
		/*Copy the inital condtion*/
		matrix<T> inital = init;
		matrix<T> J = matrix<T>::zeros(init.numel(),init.numel());
		matrix<T> fv = matrix<T>::zeros(init.numel(),1);
	
		/*Get the size of the system*/
		const int N = inital.numel();
	
		/*Set up the output file*/
		FILE *file = fopen(name,"w");
		if(!file) {std::cout << "Cannot open file\n"; exit(1);} else {}
	
		/*Counter for when to record the data*/
		unsigned int steps = (unsigned int) (end_time/dt)/(T)samples;
	
		/*Start solving*/
		unsigned int i = 0;
		while(dt*(i) <= end_time)
		{
			if(steps == (unsigned int) (end_time/dt)/(T)samples)
			{
				/*Tells how far into the program it is*/
				std::cout << 100*(dt*(double)i/(double)end_time) << "% Done\n";
			
				/*Print results to text file*/
				fprintf(file,"%.14f\t",dt*i);
				for(int i = 0; i < inital.numel(); i++)
				{
					fprintf(file,"%.14f\t",inital(i));
				}
				fprintf(file,"\n");
				steps = 0;
			}
			else {steps++;}
		
			/*Solves for the next time step*/
			/*Use newtons root finding method to solve for the future timestep*/
		
			/*Our guess values*/
			matrix<T> g1 = inital + 100.*matrix<T>::ones(inital.numel(),1);
			matrix<T> g2 = inital;
			/*Start the root finding method*/
			unsigned int count = 0;
			while(norm_2(g1 - g2)/norm_2(g1) > tol && count++ < max_iter)
			{
				/*Update values*/
				g1 = g2;
				/*The Jacobian and function vector of the resiude*/
			
				/*Initalize the Jacobian*/
				//#pragma omp parallel for
				for(int j = 0; j < N; j++)
				{
					matrix<T> copy_inital = g1;
					copy_inital(j) += step_size;
					
					J.set_col(j,-dt*(F(copy_inital,(T)dt*(i+1)) - F(g1,(T)dt*(i+1)))/step_size);
				}
				J += matrix<T>::eye(N);
			
				/*Initalize the function vector*/
				fv = g1 - inital - dt*F(g1);
				g2 = g1 - gauss<T>::solve(J,fv);
				//std::cout << count << "\n";
			}
			if(norm_2(g2 - g1)/norm_2(g2) > tol) {std::cout << "Could not find solution in desired tolerance\n";} else{}
			inital = g2;
			
			/*Move to the next time step*/
			i++;
		}
		/*Close the File and return the final solution*/
		fclose(file);
		std::cout << 100 << "% Done\n";
		return inital;
	}
	
	/*Matrix functions for speed and if the PDE is linear*/
	//static matrix<T> euler(const matrix<T>&,const matrix<T>&,matrix<T> (*)(T),T,T,unsigned int,unsigned int); 
	//static matrix<T> euler(const matrix<T>&,const matrix<T>&,matrix<T> (*)(T),T,T,char [],unsigned int,unsigned int);
	
	template<class T>
	matrix<T> ode<T>::euler(const matrix<T> &init,const matrix<T> &A,matrix<T> (*F)(T),T end_time,T dt,unsigned int samples,unsigned int band)
	{
		/*Checks to make sure the init is a vector*/
		if(init.cols() != 1) {std::cout << "Runtime Error\nFunction must be a vector\n"; exit(1);} else {}
		/*Checks to see if we can get the correct number of samples in how many steps*/
		if(samples > (unsigned int)((T)end_time/dt)) {std::cout << "Runtime Error\nNumber of samples is greater then number of runs\n";exit(1);} else {}
		/*Copy the inital condtion*/
		matrix<T> inital = init;
		
		/*Set up the output file*/
		FILE *file = fopen("output.txt","w");
		if(!file) {std::cout << "Cannot open file\n";exit(1);} else {}
		
		/*Counter for when to record the data*/
		unsigned int steps = (unsigned int) (end_time/dt)/(T)samples;
		
		/*Start solving*/
		unsigned int i = 0;
		while(dt*(i) <= end_time)
		{
			if(steps == (unsigned int) (end_time/dt)/(T)samples)
			{
				/*Tells how far into the program it is*/
				std::cout << 100*(dt*(double)i/(double)end_time) << "% Done\n";
				
				/*Print results to text file*/
				fprintf(file,"%.24f\t",dt*i);
				for(int i = 0; i < inital.numel(); i++)
				{
					fprintf(file,"%.24f\t",inital(i));
				}
				fprintf(file,"\n");
				steps = 0;
			}
			else {steps++;}
			
			/*Solves for the next time step*/
			inital += dt*(mat_mull_vec_band(A,inital,band) + F(dt*i));

			/*Move to the next time step*/
			i++;
			
		}
		/*Close the File and return the final solution*/
		fclose(file);
		std::cout << 100 << "% Done\n";
		return inital;
	}
	
	/*Implict Methods*/
	/*Must set the form in proper form where it equals a residue*/
	template<class T>
	matrix<T> ode<T>::ieuler_tridiag(const matrix<T> &init,matrix<T> (*F)(const matrix<T>&),T end_time,T dt,unsigned int samples,T tol,T step_size,unsigned int max_iter)
	{
		/*Checks to make sure the init is a vector*/
		if(init.cols() != 1) {std::cout << "Runtime Error\nFunction must be a vector\n"; exit(1);} else {}
		/*Checks to see if we can get the correct number of samples in how many steps*/
		if(samples > (unsigned int)((T)end_time/dt)) {std::cout << "Runtime Error\nNumber of samples is greater then number of runs\n";exit(1);} else {}
		/*Copy the inital condtion*/
		matrix<T> inital = init;
		matrix<T> J = matrix<T>::zeros(init.numel(),init.numel());
		matrix<T> fv = matrix<T>::zeros(init.numel(),1);
	
		/*Get the size of the system*/
		const int N = inital.numel();
	
		/*Set up the output file*/
		FILE *file = fopen("output.txt","w");
		if(!file) {std::cout << "Cannot open file\n"; exit(1);} else {}
	
		/*Counter for when to record the data*/
		unsigned int steps = (unsigned int) (end_time/dt)/(T)samples;
	
		/*Start solving*/
		unsigned int i = 0;
		while(dt*(i) <= end_time)
		{
			if(steps == (unsigned int) (end_time/dt)/(T)samples)
			{
				/*Tells how far into the program it is*/
				std::cout << 100*(dt*(double)i/(double)end_time) << "% Done\n";
			
				/*Print results to text file*/
				fprintf(file,"%.14f\t",dt*i);
				for(int i = 0; i < inital.numel(); i++)
				{
					fprintf(file,"%.14f\t",inital(i));
				}
				fprintf(file,"\n");
				steps = 0;
			}
			else {steps++;}
		
			/*Solves for the next time step*/
			/*Use newtons root finding method to solve for the future timestep*/
		
			/*Our guess values*/
			matrix<T> g1 = inital + 100.*matrix<T>::ones(inital.numel(),1);
			matrix<T> g2 = inital;
			/*Start the root finding method*/
			unsigned int count = 0;
			while(norm_2(g1 - g2)/norm_2(g1) > tol && count++ < max_iter)
			{
				/*Update values*/
				g1 = g2;
				/*The Jacobian and function vector of the resiude*/
			
				/*Initalize the Jacobian*/
				//#pragma omp parallel for
				for(int j = 0; j < N; j++)
				{
					matrix<T> copy_inital = g1;
					copy_inital(j) += step_size;
					
					J.set_col(j,-dt*(F(copy_inital) - F(g1))/step_size);
				}
				J += matrix<T>::eye(N);
			
				/*Initalize the function vector*/
				fv = g1 - inital - dt*F(g1);
				g2 = g1 - gauss<T>::solve_tridiag(J,fv);
				//std::cout << count << "\n";
			}
			if(norm_2(g2 - g1)/norm_2(g2) > tol) {std::cout << "Could not find solution in desired tolerance\n";} else{}
			inital = g2;
		
			/*Move to the next time step*/
			i++;
		}
		/*Close the File and return the final solution*/
		fclose(file);
		std::cout << 100 << "% Done\n";
		return inital;
	}
	
	template<class T>
	matrix<T> ode<T>::ieuler_tridiag(const matrix<T> &init,matrix<T> (*F)(const matrix<T>&,T),T end_time,T dt,unsigned int samples,T tol,T step_size,unsigned int max_iter)
	{
		/*Checks to make sure the init is a vector*/
		if(init.cols() != 1) {std::cout << "Runtime Error\nFunction must be a vector\n"; exit(1);} else {}
		/*Checks to see if we can get the correct number of samples in how many steps*/
		if(samples > (unsigned int)((T)end_time/dt)) {std::cout << "Runtime Error\nNumber of samples is greater then number of runs\n";exit(1);} else {}
		/*Copy the inital condtion*/
		matrix<T> inital = init;
		matrix<T> J = matrix<T>::zeros(init.numel(),init.numel());
		matrix<T> fv = matrix<T>::zeros(init.numel(),1);
	
		/*Get the size of the system*/
		const int N = inital.numel();
	
		/*Set up the output file*/
		FILE *file = fopen("output.txt","w");
		if(!file) {std::cout << "Cannot open file\n"; exit(1);} else {}
	
		/*Counter for when to record the data*/
		unsigned int steps = (unsigned int) (end_time/dt)/(T)samples;
	
		/*Start solving*/
		unsigned int i = 0;
		while(dt*(i) <= end_time)
		{
			if(steps == (unsigned int) (end_time/dt)/(T)samples)
			{
				/*Tells how far into the program it is*/
				std::cout << 100*(dt*(double)i/(double)end_time) << "% Done\n";
			
				/*Print results to text file*/
				fprintf(file,"%.14f\t",dt*i);
				for(int i = 0; i < inital.numel(); i++)
				{
					fprintf(file,"%.14f\t",inital(i));
				}
				fprintf(file,"\n");
				steps = 0;
			}
			else {steps++;}
		
			/*Solves for the next time step*/
			/*Use newtons root finding method to solve for the future timestep*/
		
			/*Our guess values*/
			matrix<T> g1 = inital + 100.*matrix<T>::ones(inital.numel(),1);
			matrix<T> g2 = inital;
			/*Start the root finding method*/
			unsigned int count = 0;
			while(norm_2(g1 - g2)/norm_2(g1) > tol && count++ < max_iter)
			{
				/*Update values*/
				g1 = g2;
				/*The Jacobian and function vector of the resiude*/
			
				/*Initalize the Jacobian*/
				//#pragma omp parallel for
				for(int j = 0; j < N; j++)
				{
					matrix<T> copy_inital = g1;
					copy_inital(j) += step_size;
					
					J.set_col(j,-dt*(F(copy_inital,(T)dt*(i+1)) - F(g1,(T)dt*(i+1)))/step_size);
				}
				J += matrix<T>::eye(N);
			
				/*Initalize the function vector*/
				fv = g1 - inital - dt*F(g1);
				g2 = g1 - gauss<T>::solve_tridiag(J,fv);
				//std::cout << count << "\n";
			}
			if(norm_2(g2 - g1)/norm_2(g2) > tol) {std::cout << "Could not find solution in desired tolerance\n";} else{}
			inital = g2;
			
			/*Move to the next time step*/
			i++;
		}
		/*Close the File and return the final solution*/
		fclose(file);
		std::cout << 100 << "% Done\n";
		return inital;
	}
	
	template<class T>
	matrix<T> ode<T>::ieuler_tridiag(const matrix<T> &init,matrix<T> (*F)(const matrix<T>&),T end_time,T dt,char name[],unsigned int samples,T tol,T step_size,unsigned int max_iter)
	{
		/*Checks to make sure the init is a vector*/
		if(init.cols() != 1) {std::cout << "Runtime Error\nFunction must be a vector\n"; exit(1);} else {}
		/*Checks to see if we can get the correct number of samples in how many steps*/
		if(samples > (unsigned int)((T)end_time/dt)) {std::cout << "Runtime Error\nNumber of samples is greater then number of runs\n";exit(1);} else {}
		/*Copy the inital condtion*/
		matrix<T> inital = init;
		matrix<T> J = matrix<T>::zeros(init.numel(),init.numel());
		matrix<T> fv = matrix<T>::zeros(init.numel(),1);
	
		/*Get the size of the system*/
		const int N = inital.numel();
	
		/*Set up the output file*/
		FILE *file = fopen(name,"w");
		if(!file) {std::cout << "Cannot open file\n"; exit(1);} else {}
	
		/*Counter for when to record the data*/
		unsigned int steps = (unsigned int) (end_time/dt)/(T)samples;
	
		/*Start solving*/
		unsigned int i = 0;
		/*Initalize the relative error*/
		T rel_error = 0.;
		while(dt*(i) <= end_time)
		{
			if(steps == (unsigned int) (end_time/dt)/(T)samples)
			{
				/*Tells how far into the program it is*/
				std::cout << 100*(dt*(double)i/(double)end_time) << "% Done\n";
			
				/*Print results to text file*/
				fprintf(file,"%.14f\t",dt*i);
				for(int i = 0; i < inital.numel(); i++)
				{
					fprintf(file,"%.14f\t",inital(i));
				}
				fprintf(file,"\n");
				steps = 0;
			}
			else {steps++;}
		
			/*Solves for the next time step*/
			/*Use newtons root finding method to solve for the future timestep*/
		
			/*Our guess values*/
			matrix<T> g1 = inital + 100.*matrix<T>::ones(inital.numel(),1);
			matrix<T> g2 = inital;
			/*Start the root finding method*/
			unsigned int count = 0;
			while(norm_2(g1 - g2)/norm_2(g1) > tol && count++ < max_iter)
			{
				/*Update values*/
				g1 = g2;
				/*The Jacobian and function vector of the resiude*/
			
				/*Initalize the Jacobian*/
				//#pragma omp parallel for
				for(int j = 0; j < N; j++)
				{
					matrix<T> copy_inital = g1;
					copy_inital(j) += step_size;
					
					J.set_col(j,-dt*(F(copy_inital) - F(g1))/step_size);
				}
				J += matrix<T>::eye(N);
			
				/*Initalize the function vector*/
				fv = g1 - inital - dt*F(g1);
				g2 = g1 - gauss<T>::solve_tridiag(J,fv);
				//std::cout << count << "\n";
			}
			if(norm_2(g2 - g1)/norm_2(g2) > tol) {std::cout << "Could not find solution in desired tolerance\n";} else{}
			inital = g2;
			
			/*Move to the next time step*/
			i++;
		}
		/*Close the File and return the final solution*/
		fclose(file);
		std::cout << 100 << "% Done\n";
		return inital;
	}
	
	template<class T>
	matrix<T> ode<T>::ieuler_tridiag(const matrix<T> &init,matrix<T> (*F)(const matrix<T>&,T),T end_time,T dt,char name[],unsigned int samples,T tol,T step_size,unsigned int max_iter)
	{
		/*Checks to make sure the init is a vector*/
		if(init.cols() != 1) {std::cout << "Runtime Error\nFunction must be a vector\n"; exit(1);} else {}
		/*Checks to see if we can get the correct number of samples in how many steps*/
		if(samples > (unsigned int)((T)end_time/dt)) {std::cout << "Runtime Error\nNumber of samples is greater then number of runs\n";exit(1);} else {}
		/*Copy the inital condtion*/
		matrix<T> inital = init;
		matrix<T> J = matrix<T>::zeros(init.numel(),init.numel());
		matrix<T> fv = matrix<T>::zeros(init.numel(),1);
	
		/*Get the size of the system*/
		const int N = inital.numel();
	
		/*Set up the output file*/
		FILE *file = fopen(name,"w");
		if(!file) {std::cout << "Cannot open file\n"; exit(1);} else {}
	
		/*Counter for when to record the data*/
		unsigned int steps = (unsigned int) (end_time/dt)/(T)samples;
	
		/*Start solving*/
		unsigned int i = 0;
		while(dt*(i) <= end_time)
		{
			if(steps == (unsigned int) (end_time/dt)/(T)samples)
			{
				/*Tells how far into the program it is*/
				std::cout << 100*(dt*(double)i/(double)end_time) << "% Done\n";
			
				/*Print results to text file*/
				fprintf(file,"%.14f\t",dt*i);
				for(int i = 0; i < inital.numel(); i++)
				{
					fprintf(file,"%.14f\t",inital(i));
				}
				fprintf(file,"\n");
				steps = 0;
			}
			else {steps++;}
		
			/*Solves for the next time step*/
			/*Use newtons root finding method to solve for the future timestep*/
		
			/*Our guess values*/
			matrix<T> g1 = inital + 100.*matrix<T>::ones(inital.numel(),1);
			matrix<T> g2 = inital;
			/*Start the root finding method*/
			unsigned int count = 0;
			while(norm_2(g1 - g2)/norm_2(g1) > tol && count++ < max_iter)
			{
				/*Update values*/
				g1 = g2;
				/*The Jacobian and function vector of the resiude*/
			
				/*Initalize the Jacobian*/
				//#pragma omp parallel for
				for(int j = 0; j < N; j++)
				{
					matrix<T> copy_inital = g1;
					copy_inital(j) += step_size;
					
					J.set_col(j,-dt*(F(copy_inital,(T)dt*(i+1)) - F(g1,(T)dt*(i+1)))/step_size);
				}
				J += matrix<T>::eye(N);
			
				/*Initalize the function vector*/
				fv = g1 - inital - dt*F(g1);
				g2 = g1 - gauss<T>::solve_tridiag(J,fv);
				//std::cout << count << "\n";
			}
			if(norm_2(g2 - g1)/norm_2(g2) > tol) {std::cout << "Could not find solution in desired tolerance\n";} else{}
			inital = g2;
			
			/*Move to the next time step*/
			i++;
		}
		/*Close the File and return the final solution*/
		fclose(file);
		std::cout << 100 << "% Done\n";
		return inital;
	}
	
	/*Implict Methods*/
	/*Must set the form in proper form where it equals a residue*/
	template<class T>
	matrix<T> ode<T>::crank_nicolson(const matrix<T> &init,matrix<T> (*F)(const matrix<T>&),T end_time,T dt,unsigned int samples,T tol,T step_size,unsigned int max_iter)
	{
		/*Checks to make sure the init is a vector*/
		if(init.cols() != 1) {std::cout << "Runtime Error\nFunction must be a vector\n"; exit(1);} else {}
		/*Checks to see if we can get the correct number of samples in how many steps*/
		if(samples > (unsigned int)((T)end_time/dt)) {std::cout << "Runtime Error\nNumber of samples is greater then number of runs\n";exit(1);} else {}
		/*Copy the inital condtion*/
		matrix<T> inital = init;
		matrix<T> J = matrix<T>::zeros(init.numel(),init.numel());
		matrix<T> fv = matrix<T>::zeros(init.numel(),1);
	
		/*Get the size of the system*/
		const int N = inital.numel();
	
		/*Set up the output file*/
		FILE *file = fopen("output.txt","w");
		if(!file) {std::cout << "Cannot open file\n"; exit(1);} else {}
	
		/*Counter for when to record the data*/
		unsigned int steps = (unsigned int) (end_time/dt)/(T)samples;
	
		/*Start solving*/
		unsigned int i = 0;
		while(dt*(i) <= end_time)
		{
			if(steps == (unsigned int) (end_time/dt)/(T)samples)
			{
				/*Tells how far into the program it is*/
				std::cout << 100*(dt*(double)i/(double)end_time) << "% Done\n";
			
				/*Print results to text file*/
				fprintf(file,"%.14f\t",dt*i);
				for(int i = 0; i < inital.numel(); i++)
				{
					fprintf(file,"%.14f\t",inital(i));
				}
				fprintf(file,"\n");
				steps = 0;
			}
			else {steps++;}
		
			/*Solves for the next time step*/
			/*Use newtons root finding method to solve for the future timestep*/
		
			/*Our guess values*/
			matrix<T> g1 = inital + 100.*matrix<T>::ones(inital.numel(),1);
			matrix<T> g2 = inital;
			/*Start the root finding method*/
			unsigned int count = 0;
			while(norm_2(g1 - g2)/norm_2(g1) > tol && count++ < max_iter)
			{
				/*Update values*/
				g1 = g2;
				/*The Jacobian and function vector of the resiude*/
			
				/*Initalize the Jacobian*/
				//#pragma omp parallel for
				for(int j = 0; j < N; j++)
				{
					matrix<T> copy_inital = g1;
					copy_inital(j) += step_size;
					
					J.set_col(j,-.5*dt*(F(copy_inital) - F(g1))/step_size);
				}
				J += matrix<T>::eye(N);
			
				/*Initalize the function vector*/
				fv = g1 - inital - .5*dt*(F(g1) + F(inital));
				g2 = g1 - gauss<T>::solve(J,fv);
				//std::cout << count << "\n";
			}
			if(norm_2(g2 - g1)/norm_2(g2) > tol) {std::cout << "Could not find solution in desired tolerance\n";} else{}
			inital = g2;
			/*Move to next time step*/
			i++;
		}
		/*Close the File and return the final solution*/
		fclose(file);
		std::cout << 100 << "% Done\n";
		return inital;
	}
	
	template<class T>
	matrix<T> ode<T>::crank_nicolson(const matrix<T> &init,matrix<T> (*F)(const matrix<T>&,T),T end_time,T dt,unsigned int samples,T tol,T step_size,unsigned int max_iter)
	{
		/*Checks to make sure the init is a vector*/
		if(init.cols() != 1) {std::cout << "Runtime Error\nFunction must be a vector\n"; exit(1);} else {}
		/*Checks to see if we can get the correct number of samples in how many steps*/
		if(samples > (unsigned int)((T)end_time/dt)) {std::cout << "Runtime Error\nNumber of samples is greater then number of runs\n";exit(1);} else {}
		/*Copy the inital condtion*/
		matrix<T> inital = init;
		matrix<T> J = matrix<T>::zeros(init.numel(),init.numel());
		matrix<T> fv = matrix<T>::zeros(init.numel(),1);
	
		/*Get the size of the system*/
		const int N = inital.numel();
	
		/*Set up the output file*/
		FILE *file = fopen("output.txt","w");
		if(!file) {std::cout << "Cannot open file\n"; exit(1);} else {}
	
		/*Counter for when to record the data*/
		unsigned int steps = (unsigned int) (end_time/dt)/(T)samples;
	
		/*Start solving*/
		unsigned int i = 0;
		while(dt*(i) <= end_time)
		{
			if(steps == (unsigned int) (end_time/dt)/(T)samples)
			{
				/*Tells how far into the program it is*/
				std::cout << 100*(dt*(double)i/(double)end_time) << "% Done\n";
			
				/*Print results to text file*/
				fprintf(file,"%.14f\t",dt*i);
				for(int i = 0; i < inital.numel(); i++)
				{
					fprintf(file,"%.14f\t",inital(i));
				}
				fprintf(file,"\n");
				steps = 0;
			}
			else {steps++;}
		
			/*Solves for the next time step*/
			/*Use newtons root finding method to solve for the future timestep*/
		
			/*Our guess values*/
			matrix<T> g1 = inital + 100.*matrix<T>::ones(inital.numel(),1);
			matrix<T> g2 = inital;
			/*Start the root finding method*/
			unsigned int count = 0;
			while(norm_2(g1 - g2)/norm_2(g1) > tol && count++ < max_iter)
			{
				/*Update values*/
				g1 = g2;
				/*The Jacobian and function vector of the resiude*/
			
				/*Initalize the Jacobian*/
				//#pragma omp parallel for
				for(int j = 0; j < N; j++)
				{
					matrix<T> copy_inital = g1;
					copy_inital(j) += step_size;
					
					J.set_col(j,-.5*dt*(F(copy_inital,(T)dt*(i+1)) - F(g1,(T)dt*(i+1)))/step_size);
				}
				J += matrix<T>::eye(N);
			
				/*Initalize the function vector*/
				fv = g1 - inital - .5*dt*(F(g1,dt*(i+1)) + F(inital,dt*i));
				g2 = g1 - gauss<T>::solve(J,fv);
				//std::cout << count << "\n";
			}
			if(norm_2(g2 - g1)/norm_2(g2) > tol) {std::cout << "Could not find solution in desired tolerance\n";} else{}
			inital = g2;
			
			/*Move to the next time step*/
			i++;
		}
		/*Close the File and return the final solution*/
		fclose(file);
		std::cout << 100 << "% Done\n";
		return inital;
	}
	
	template<class T>
	matrix<T> ode<T>::crank_nicolson(const matrix<T> &init,matrix<T> (*F)(const matrix<T>&),T end_time,T dt,char name[],unsigned int samples,T tol,T step_size,unsigned int max_iter)
	{
		/*Checks to make sure the init is a vector*/
		if(init.cols() != 1) {std::cout << "Runtime Error\nFunction must be a vector\n"; exit(1);} else {}
		/*Checks to see if we can get the correct number of samples in how many steps*/
		if(samples > (unsigned int)((T)end_time/dt)) {std::cout << "Runtime Error\nNumber of samples is greater then number of runs\n";exit(1);} else {}
		/*Copy the inital condtion*/
		matrix<T> inital = init;
		matrix<T> J = matrix<T>::zeros(init.numel(),init.numel());
		matrix<T> fv = matrix<T>::zeros(init.numel(),1);
	
		/*Get the size of the system*/
		const int N = inital.numel();
	
		/*Set up the output file*/
		FILE *file = fopen(name,"w");
		if(!file) {std::cout << "Cannot open file\n"; exit(1);} else {}
	
		/*Counter for when to record the data*/
		unsigned int steps = (unsigned int) (end_time/dt)/(T)samples;
	
		/*Start solving*/
		unsigned int i = 0;
		while(dt*(i) <= end_time)
		{
			if(steps == (unsigned int) (end_time/dt)/(T)samples)
			{
				/*Tells how far into the program it is*/
				std::cout << 100*(dt*(double)i/(double)end_time) << "% Done\n";
			
				/*Print results to text file*/
				fprintf(file,"%.14f\t",dt*i);
				for(int i = 0; i < inital.numel(); i++)
				{
					fprintf(file,"%.14f\t",inital(i));
				}
				fprintf(file,"\n");
				steps = 0;
			}
			else {steps++;}
		
			/*Solves for the next time step*/
			/*Use newtons root finding method to solve for the future timestep*/
		
			/*Our guess values*/
			matrix<T> g1 = inital + 100.*matrix<T>::ones(inital.numel(),1);
			matrix<T> g2 = inital;
			/*Start the root finding method*/
			unsigned int count = 0;
			while(norm_2(g1 - g2)/norm_2(g1) > tol && count++ < max_iter)
			{
				/*Update values*/
				g1 = g2;
				/*The Jacobian and function vector of the resiude*/
			
				/*Initalize the Jacobian*/
				//#pragma omp parallel for
				for(int j = 0; j < N; j++)
				{
					matrix<T> copy_inital = g1;
					copy_inital(j) += step_size;
					
					J.set_col(j,-.5*dt*(F(copy_inital) - F(g1))/step_size);
				}
				J += matrix<T>::eye(N);
			
				/*Initalize the function vector*/
				fv = g1 - inital - .5*dt*(F(g1)+F(inital));
				g2 = g1 - gauss<T>::solve(J,fv);
				//std::cout << count << "\n";
			}
			if(norm_2(g2 - g1)/norm_2(g2) > tol) {std::cout << "Could not find solution in desired tolerance\n";} else{}
			inital = g2;
			/*Move to the next time step*/
			i++;
		}
		/*Close the File and return the final solution*/
		fclose(file);
		std::cout << 100 << "% Done\n";
		return inital;
	}
	
	template<class T>
	matrix<T> ode<T>::crank_nicolson(const matrix<T> &init,matrix<T> (*F)(const matrix<T>&,T),T end_time,T dt,char name[],unsigned int samples,T tol,T step_size,unsigned int max_iter)
	{
		/*Checks to make sure the init is a vector*/
		if(init.cols() != 1) {std::cout << "Runtime Error\nFunction must be a vector\n"; exit(1);} else {}
		/*Checks to see if we can get the correct number of samples in how many steps*/
		if(samples > (unsigned int)((T)end_time/dt)) {std::cout << "Runtime Error\nNumber of samples is greater then number of runs\n";exit(1);} else {}
		/*Copy the inital condtion*/
		matrix<T> inital = init;
		matrix<T> J = matrix<T>::zeros(init.numel(),init.numel());
		matrix<T> fv = matrix<T>::zeros(init.numel(),1);
	
		/*Get the size of the system*/
		const int N = inital.numel();
	
		/*Set up the output file*/
		FILE *file = fopen(name,"w");
		if(!file) {std::cout << "Cannot open file\n"; exit(1);} else {}
	
		/*Counter for when to record the data*/
		unsigned int steps = (unsigned int) (end_time/dt)/(T)samples;
	
		/*Start solving*/
		unsigned int i = 0;
		while(dt*(i) <= end_time)
		{
			if(steps == (unsigned int) (end_time/dt)/(T)samples)
			{
				/*Tells how far into the program it is*/
				std::cout << 100*(dt*(double)i/(double)end_time) << "% Done\n";
			
				/*Print results to text file*/
				fprintf(file,"%.14f\t",dt*i);
				for(int i = 0; i < inital.numel(); i++)
				{
					fprintf(file,"%.14f\t",inital(i));
				}
				fprintf(file,"\n");
				steps = 0;
			}
			else {steps++;}
		
			/*Solves for the next time step*/
			/*Use newtons root finding method to solve for the future timestep*/
		
			/*Our guess values*/
			matrix<T> g1 = inital + 100.*matrix<T>::ones(inital.numel(),1);
			matrix<T> g2 = inital;
			/*Start the root finding method*/
			unsigned int count = 0;
			while(norm_2(g1 - g2)/norm_2(g1) > tol && count++ < max_iter)
			{
				/*Update values*/
				g1 = g2;
				/*The Jacobian and function vector of the resiude*/
			
				/*Initalize the Jacobian*/
				//#pragma omp parallel for
				for(int j = 0; j < N; j++)
				{
					matrix<T> copy_inital = g1;
					copy_inital(j) += step_size;
					
					J.set_col(j,-.5*dt*(F(copy_inital,(T)dt*(i+1)) - F(g1,(T)dt*(i+1)))/step_size);
				}
				J += matrix<T>::eye(N);
			
				/*Initalize the function vector*/
				fv = g1 - inital - .5*dt*(F(g1)+F(inital));
				g2 = g1 - gauss<T>::solve(J,fv);
				//std::cout << count << "\n";
			}
			if(norm_2(g2 - g1)/norm_2(g2) > tol) {std::cout << "Could not find solution in desired tolerance\n";} else{}
			inital = g2;
			/*Move to the next time step*/
			i++;
		}
		/*Close the File and return the final solution*/
		fclose(file);
		std::cout << 100 << "% Done\n";
		return inital;
	}
	
	template<class T>
	matrix<T> ode<T>::crank_nicolson_tridiag(const matrix<T> &init,matrix<T> (*F)(const matrix<T>&),T end_time,T dt,unsigned int samples,T tol,T step_size,unsigned int max_iter)
	{
		/*Checks to make sure the init is a vector*/
		if(init.cols() != 1) {std::cout << "Runtime Error\nFunction must be a vector\n"; exit(1);} else {}
		/*Checks to see if we can get the correct number of samples in how many steps*/
		if(samples > (unsigned int)((T)end_time/dt)) {std::cout << "Runtime Error\nNumber of samples is greater then number of runs\n";exit(1);} else {}
		/*Copy the inital condtion*/
		matrix<T> inital = init;
		matrix<T> J = matrix<T>::zeros(init.numel(),init.numel());
		matrix<T> fv = matrix<T>::zeros(init.numel(),1);
	
		/*Get the size of the system*/
		const int N = inital.numel();
	
		/*Set up the output file*/
		FILE *file = fopen("output.txt","w");
		if(!file) {std::cout << "Cannot open file\n"; exit(1);} else {}
	
		/*Counter for when to record the data*/
		unsigned int steps = (unsigned int) (end_time/dt)/(T)samples;
	
		/*Start solving*/
		unsigned int i = 0;
		while(dt*(i) <= end_time)
		{
			if(steps == (unsigned int) (end_time/dt)/(T)samples)
			{
				/*Tells how far into the program it is*/
				std::cout << 100*(dt*(double)i/(double)end_time) << "% Done\n";
			
				/*Print results to text file*/
				fprintf(file,"%.14f\t",dt*i);
				for(int i = 0; i < inital.numel(); i++)
				{
					fprintf(file,"%.14f\t",inital(i));
				}
				fprintf(file,"\n");
				steps = 0;
			}
			else {steps++;}
		
			/*Solves for the next time step*/
			/*Use newtons root finding method to solve for the future timestep*/
		
			/*Our guess values*/
			matrix<T> g1 = inital + 100.*matrix<T>::ones(inital.numel(),1);
			matrix<T> g2 = inital;
			/*Start the root finding method*/
			unsigned int count = 0;
			while(norm_2(g1 - g2)/norm_2(g1) > tol && count++ < max_iter)
			{
				/*Update values*/
				g1 = g2;
				/*The Jacobian and function vector of the resiude*/
			
				/*Initalize the Jacobian*/
				//#pragma omp parallel for
				for(int j = 0; j < N; j++)
				{
					matrix<T> copy_inital = g1;
					copy_inital(j) += step_size;
					
					J.set_col(j,-.5*dt*(F(copy_inital) - F(g1))/step_size);
				}
				J += matrix<T>::eye(N);
			
				/*Initalize the function vector*/
				fv = g1 - inital - .5*dt*(F(g1) + F(inital));
				g2 = g1 - gauss<T>::solve_tridiag(J,fv);
				//std::cout << count << "\n";
			}
			if(norm_2(g2 - g1)/norm_2(g2) > tol) {std::cout << "Could not find solution in desired tolerance\n";} else{}
			inital = g2;
			/*Move to next time step*/
			i++;
		}
		/*Close the File and return the final solution*/
		fclose(file);
		std::cout << 100 << "% Done\n";
		return inital;
	}
	
	template<class T>
	matrix<T> ode<T>::crank_nicolson_tridiag(const matrix<T> &init,matrix<T> (*F)(const matrix<T>&,T),T end_time,T dt,unsigned int samples,T tol,T step_size,unsigned int max_iter)
	{
		/*Checks to make sure the init is a vector*/
		if(init.cols() != 1) {std::cout << "Runtime Error\nFunction must be a vector\n"; exit(1);} else {}
		/*Checks to see if we can get the correct number of samples in how many steps*/
		if(samples > (unsigned int)((T)end_time/dt)) {std::cout << "Runtime Error\nNumber of samples is greater then number of runs\n";exit(1);} else {}
		/*Copy the inital condtion*/
		matrix<T> inital = init;
		matrix<T> J = matrix<T>::zeros(init.numel(),init.numel());
		matrix<T> fv = matrix<T>::zeros(init.numel(),1);
	
		/*Get the size of the system*/
		const int N = inital.numel();
	
		/*Set up the output file*/
		FILE *file = fopen("output.txt","w");
		if(!file) {std::cout << "Cannot open file\n"; exit(1);} else {}
	
		/*Counter for when to record the data*/
		unsigned int steps = (unsigned int) (end_time/dt)/(T)samples;
	
		/*Start solving*/
		unsigned int i = 0;
		while(dt*(i) <= end_time)
		{
			if(steps == (unsigned int) (end_time/dt)/(T)samples)
			{
				/*Tells how far into the program it is*/
				std::cout << 100*(dt*(double)i/(double)end_time) << "% Done\n";
			
				/*Print results to text file*/
				fprintf(file,"%.14f\t",dt*i);
				for(int i = 0; i < inital.numel(); i++)
				{
					fprintf(file,"%.14f\t",inital(i));
				}
				fprintf(file,"\n");
				steps = 0;
			}
			else {steps++;}
		
			/*Solves for the next time step*/
			/*Use newtons root finding method to solve for the future timestep*/
		
			/*Our guess values*/
			matrix<T> g1 = inital + 100.*matrix<T>::ones(inital.numel(),1);
			matrix<T> g2 = inital;
			/*Start the root finding method*/
			unsigned int count = 0;
			while(norm_2(g1 - g2)/norm_2(g1) > tol && count++ < max_iter)
			{
				/*Update values*/
				g1 = g2;
				/*The Jacobian and function vector of the resiude*/
			
				/*Initalize the Jacobian*/
				//#pragma omp parallel for
				for(int j = 0; j < N; j++)
				{
					matrix<T> copy_inital = g1;
					copy_inital(j) += step_size;
					
					J.set_col(j,-.5*dt*(F(copy_inital,(T)dt*(i+1)) - F(g1,(T)dt*(i+1)))/step_size);
				}
				J += matrix<T>::eye(N);
			
				/*Initalize the function vector*/
				fv = g1 - inital - .5*dt*(F(g1,dt*(i+1)) + F(inital,dt*i));
				g2 = g1 - gauss<T>::solve_tridiag(J,fv);
				//std::cout << count << "\n";
			}
			if(norm_2(g2 - g1)/norm_2(g2) > tol) {std::cout << "Could not find solution in desired tolerance\n";} else{}
			inital = g2;
			
			/*Move to the next time step*/
			i++;
		}
		/*Close the File and return the final solution*/
		fclose(file);
		std::cout << 100 << "% Done\n";
		return inital;
	}
	
	template<class T>
	matrix<T> ode<T>::crank_nicolson_tridiag(const matrix<T> &init,matrix<T> (*F)(const matrix<T>&),T end_time,T dt,char name[],unsigned int samples,T tol,T step_size,unsigned int max_iter)
	{
		/*Checks to make sure the init is a vector*/
		if(init.cols() != 1) {std::cout << "Runtime Error\nFunction must be a vector\n"; exit(1);} else {}
		/*Checks to see if we can get the correct number of samples in how many steps*/
		if(samples > (unsigned int)((T)end_time/dt)) {std::cout << "Runtime Error\nNumber of samples is greater then number of runs\n";exit(1);} else {}
		/*Copy the inital condtion*/
		matrix<T> inital = init;
		matrix<T> J = matrix<T>::zeros(init.numel(),init.numel());
		matrix<T> fv = matrix<T>::zeros(init.numel(),1);
	
		/*Get the size of the system*/
		const int N = inital.numel();
	
		/*Set up the output file*/
		FILE *file = fopen(name,"w");
		if(!file) {std::cout << "Cannot open file\n"; exit(1);} else {}
	
		/*Counter for when to record the data*/
		unsigned int steps = (unsigned int) (end_time/dt)/(T)samples;
	
		/*Start solving*/
		unsigned int i = 0;
		while(dt*(i) <= end_time)
		{
			if(steps == (unsigned int) (end_time/dt)/(T)samples)
			{
				/*Tells how far into the program it is*/
				std::cout << 100*(dt*(double)i/(double)end_time) << "% Done\n";
			
				/*Print results to text file*/
				fprintf(file,"%.14f\t",dt*i);
				for(int i = 0; i < inital.numel(); i++)
				{
					fprintf(file,"%.14f\t",inital(i));
				}
				fprintf(file,"\n");
				steps = 0;
			}
			else {steps++;}
		
			/*Solves for the next time step*/
			/*Use newtons root finding method to solve for the future timestep*/
		
			/*Our guess values*/
			matrix<T> g1 = inital + 100.*matrix<T>::ones(inital.numel(),1);
			matrix<T> g2 = inital;
			/*Start the root finding method*/
			unsigned int count = 0;
			while(norm_2(g1 - g2)/norm_2(g1) > tol && count++ < max_iter)
			{
				/*Update values*/
				g1 = g2;
				/*The Jacobian and function vector of the resiude*/
			
				/*Initalize the Jacobian*/
				//#pragma omp parallel for
				for(int j = 0; j < N; j++)
				{
					matrix<T> copy_inital = g1;
					copy_inital(j) += step_size;
					
					J.set_col(j,-.5*dt*(F(copy_inital) - F(g1))/step_size);
				}
				J += matrix<T>::eye(N);
			
				/*Initalize the function vector*/
				fv = g1 - inital - .5*dt*(F(g1)+F(inital));
				g2 = g1 - gauss<T>::solve_tridiag(J,fv);
				//std::cout << count << "\n";
			}
			if(norm_2(g2 - g1)/norm_2(g2) > tol) {std::cout << "Could not find solution in desired tolerance\n";} else{}
			inital = g2;
			/*Move to the next time step*/
			i++;
		}
		/*Close the File and return the final solution*/
		fclose(file);
		std::cout << 100 << "% Done\n";
		return inital;
	}
	
	template<class T>
	matrix<T> ode<T>::crank_nicolson_tridiag(const matrix<T> &init,matrix<T> (*F)(const matrix<T>&,T),T end_time,T dt,char name[],unsigned int samples,T tol,T step_size,unsigned int max_iter)
	{
		/*Checks to make sure the init is a vector*/
		if(init.cols() != 1) {std::cout << "Runtime Error\nFunction must be a vector\n"; exit(1);} else {}
		/*Checks to see if we can get the correct number of samples in how many steps*/
		if(samples > (unsigned int)((T)end_time/dt)) {std::cout << "Runtime Error\nNumber of samples is greater then number of runs\n";exit(1);} else {}
		/*Copy the inital condtion*/
		matrix<T> inital = init;
		matrix<T> J = matrix<T>::zeros(init.numel(),init.numel());
		matrix<T> fv = matrix<T>::zeros(init.numel(),1);
	
		/*Get the size of the system*/
		const int N = inital.numel();
	
		/*Set up the output file*/
		FILE *file = fopen(name,"w");
		if(!file) {std::cout << "Cannot open file\n"; exit(1);} else {}
	
		/*Counter for when to record the data*/
		unsigned int steps = (unsigned int) (end_time/dt)/(T)samples;
	
		/*Start solving*/
		unsigned int i = 0;
		while(dt*(i) <= end_time)
		{
			if(steps == (unsigned int) (end_time/dt)/(T)samples)
			{
				/*Tells how far into the program it is*/
				std::cout << 100*(dt*(double)i/(double)end_time) << "% Done\n";
			
				/*Print results to text file*/
				fprintf(file,"%.14f\t",dt*i);
				for(int i = 0; i < inital.numel(); i++)
				{
					fprintf(file,"%.14f\t",inital(i));
				}
				fprintf(file,"\n");
				steps = 0;
			}
			else {steps++;}
		
			/*Solves for the next time step*/
			/*Use newtons root finding method to solve for the future timestep*/
		
			/*Our guess values*/
			matrix<T> g1 = inital + 100.*matrix<T>::ones(inital.numel(),1);
			matrix<T> g2 = inital;
			/*Start the root finding method*/
			unsigned int count = 0;
			while(norm_2(g1 - g2)/norm_2(g1) > tol && count++ < max_iter)
			{
				/*Update values*/
				g1 = g2;
				/*The Jacobian and function vector of the resiude*/
			
				/*Initalize the Jacobian*/
				//#pragma omp parallel for
				for(int j = 0; j < N; j++)
				{
					matrix<T> copy_inital = g1;
					copy_inital(j) += step_size;
					
					J.set_col(j,-.5*dt*(F(copy_inital,(T)dt*(i+1)) - F(g1,(T)dt*(i+1)))/step_size);
				}
				J += matrix<T>::eye(N);
			
				/*Initalize the function vector*/
				fv = g1 - inital - .5*dt*(F(g1)+F(inital));
				g2 = g1 - gauss<T>::solve_tridiag(J,fv);
				//std::cout << count << "\n";
			}
			if(norm_2(g2 - g1)/norm_2(g2) > tol) {std::cout << "Could not find solution in desired tolerance\n";} else{}
			inital = g2;
			/*Move to the next time step*/
			i++;
		}
		/*Close the File and return the final solution*/
		fclose(file);
		std::cout << 100 << "% Done\n";
		return inital;
	}
}

#define dode ode<double>
#define fode ode<float>
#define iode ode<int>

