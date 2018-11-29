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
	};

	/*Euler*/
	template<class T>
	matrix<T> ode<T>::euler(const matrix<T> &init,matrix<T> (*F)(const matrix<T>&,T),T end_time,T dt,unsigned int samples) //Techincally starts at t = 0
	{
		/*Checks to make sure the init is a vector*/
		if(init.cols() != 1) {std::cout << "Runtime Error\nFunction must be a vector\n"; exit(1);} else {}
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
			inital += dt*F(inital,dt*i);
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
			inital += dt*F(inital);
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
			inital += dt*F(inital,dt*i);
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
			inital += dt*F(inital);
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
			matrix<T> k1 = inital + .5*dt*F(inital,dt*i);
			inital += dt*F(k1,dt*i + .5*dt);
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
			matrix<T> k1 = inital + .5*dt*F(inital);
			inital += dt*F(k1);
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
			matrix<T> k1 = inital + .5*dt*F(inital,dt*i);
			inital += dt*F(k1,dt*i + .5*dt);
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
			matrix<T> k1 = inital + .5*dt*F(inital);
			inital += dt*F(k1);
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
			matrix<T> k1 = dt*F(inital);
			matrix<T> k2 = dt*F(inital + k1/(T)2);
			matrix<T> k3 = dt*F(inital - k1 + (T)2*k2);
			inital += (T)(1./6.)*(k1 + (T)4*k2 + k3);
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
			matrix<T> k1 = dt*F(inital,dt*i);
			matrix<T> k2 = dt*F(inital + k1/(T)2,dt*i + dt/(T)2);
			matrix<T> k3 = dt*F(inital - k1 + (T)2*k2,dt*i + dt);
			inital += (T)(1./6.)*(k1 + (T)4*k2 + k3);
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
			matrix<T> k1 = dt*F(inital);
			matrix<T> k2 = dt*F(inital + k1/(T)2);
			matrix<T> k3 = dt*F(inital - k1 + (T)2*k2);
			inital += (T)(1./6.)*(k1 + (T)4*k2 + k3);
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
			matrix<T> k1 = dt*F(inital,dt*i);
			matrix<T> k2 = dt*F(inital + k1/(T)2,dt*i + dt/(T)2);
			matrix<T> k3 = dt*F(inital - k1 + (T)2*k2,dt*i + dt);
			inital += (T)(1./6.)*(k1 + (T)4*k2 + k3);
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
			k1 = F(inital);
			k2 = F(inital + (dt/(T)2)*k1);
			k3 = F(inital + (dt/(T)2)*k2);
			k4 = F(inital + dt*k3);
			inital += dt*(T)(1./6.)*(k1 + (T)2*k2 + (T)2*k3 + k4);
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
			k1 = F(inital,dt*i);
			k2 = F(inital + (dt/(T)2)*k1,dt*i + dt/(T)2);
			k3 = F(inital + (dt/(T)2)*k2,dt*i + dt/(T)2);
			k4 = F(inital + dt*k3,dt*i + dt);
			inital += dt*(T)(1./6.)*(k1 + (T)2*k2 + (T)2*k3 + k4);
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
			k1 = F(inital);
			k2 = F(inital + (dt/(T)2)*k1);
			k3 = F(inital + (dt/(T)2)*k2);
			k4 = F(inital + dt*k3);
			inital += dt*(T)(1./6.)*(k1 + (T)2*k2 + (T)2*k3 + k4);
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
			k1 = F(inital,dt*i);
			k2 = F(inital + (dt/(T)2)*k1,dt*i + dt/(T)2);
			k3 = F(inital + (dt/(T)2)*k2,dt*i + dt/(T)2);
			k4 = F(inital + dt*k3,dt*i + dt);
			inital += dt*(T)(1./6.)*(k1 + (T)2*k2 + (T)2*k3 + k4);
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

