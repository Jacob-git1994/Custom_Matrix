#pragma once
#include "matrix.h"
#include "gauss.h"
#include "jacobi.h"
#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>

namespace nicole
{
	template<class T>
	class laplace
	{
	private:
		/*Set up matrix and grid*/
		matrix<T> MATRIX,FUNCTION_VECTOR,GRID,BOUNDARY_CORRECTION;
		T h;
		
		/*Forcing and boundary functions*/
		T (*Forcing)(T,T);
		matrix<T (*)(T)> boundary;
		
		/*Types of boundaries*/
		matrix<char> boundary_type;
	public:
		/*Constructors*/
		laplace(void);
		laplace(unsigned int);
		laplace(unsigned int,T,T);
		~laplace(void);
		/*Set up functions*/
		void set_forcing(T (*)(T,T));
		void set_boundary(T (*)(T),T (*)(T),T (*)(T),T (*)(T));
		void set_boundary_type(char,char,char,char);
		/*Set up the grid*/
		void set_grid(T,T);
		void set_grid_size(unsigned int);
		/*Set up problem*/
		void initalize(void);
		/*Solve*/
		matrix<T> lu_solve(void);
		matrix<T> j_solve(unsigned int,T);
	};
	
	template<class T>
	laplace<T>::laplace(void)
	{
		/*Initalize the boundary functions*/
		boundary = matrix<T (*)(T)>::shape(1,4);
		boundary_type = matrix<char>::shape(1,4);
		Forcing = nullptr;
		h = 0;
		std::cout << "Laplace Object Initalized\n";
	}
	
	template<class T>
	laplace<T>::laplace(unsigned int N)
	{
		/*Initalize the boundary functions*/
		boundary = matrix<T (*)(T)>::shape(1,4);
		boundary_type = matrix<char>::shape(1,4);
		Forcing = nullptr;
		
		/*Initalize our matrix systems*/
		MATRIX = matrix<T>::zeros(N*N,N*N);
		FUNCTION_VECTOR = matrix<T>::zeros(N*N,1);
		BOUNDARY_CORRECTION = matrix<T>::zeros(N*N,1);
		GRID = matrix<T>::linspace((T)0.,(T)1.,N+2);
		h = std::fabs(GRID(3)-GRID(2));
		std::cout << "Initalized\n";
	}

	template<class T>
	laplace<T>::laplace(unsigned int N,T left,T right)
	{
		/*Initalize the boundary functions*/
		boundary = matrix<T (*)(T)>::shape(1,4);
		boundary_type = matrix<char>::shape(1,4);
		Forcing = nullptr;
		
		/*Initalize our matrix systems*/
		MATRIX = matrix<T>::zeros(N*N,N*N);
		FUNCTION_VECTOR = matrix<T>::zeros(N*N,1);
		BOUNDARY_CORRECTION = matrix<T>::zeros(N*N,1);
		GRID = matrix<T>::linspace(left,right,N+2);
		h = std::fabs(GRID(3)-GRID(2));
		std::cout << "Initalized\n";
	}
	
	template<class T>
	laplace<T>::~laplace(void)
	{
		/*Destructors apart of the matrix class*/
		std::cout << "Laplace Solver Finished\n";
	}
	
	/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	template<class T>
	void laplace<T>::set_forcing(T (*force)(T,T))
	{
		Forcing = force;
		std::cout << "Forcing Function Initalized\n";
	}

	template<class T>
	void laplace<T>::set_boundary(T (*l)(T),T (*r)(T),T (*t)(T),T (*b)(T))
	{
		boundary(0) = l;
		boundary(1) = r;
		boundary(2) = t;
		boundary(3) = b;
		std::cout << "Boundary Functions Initalized\n";
		
	}

	template<class T>
	void laplace<T>::set_boundary_type(char l,char r,char t,char b)
	{
		boundary_type(0) = l;
		boundary_type(1) = r;
		boundary_type(2) = t;
		boundary_type(3) = b;
		std::cout << "Boundary Types Initalized\n";
	}

	template<class T>
	void laplace<T>::set_grid(T left,T right)
	{
		GRID = matrix<T>::linspace(left,right,GRID.numel());
		h = std::fabs(GRID(3)- GRID(2));
		std::cout << "Grid Rescaled\n";
	}

	template<class T>
	void laplace<T>::set_grid_size(unsigned int N)
	{
		MATRIX = matrix<T>::zeros(N*N,N*N);
		FUNCTION_VECTOR = matrix<T>::zeros(N*N,1);
		BOUNDARY_CORRECTION = matrix<T>::zeros(N*N,1);
		GRID = matrix<T>::linspace(GRID(0),GRID(GRID.numel()-1),N+2);
		h = std::fabs(GRID(3) - GRID(2));
		std::cout << "Grid Rescaled\n";
	}

	/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	/*Initalze everything and the system*/
	template<class T>
	void laplace<T>::initalize(void)
	{
		const unsigned int N = GRID.numel()-2;
		/*Initalize the forcing function*/
		unsigned int pos = 0;
#pragma omp for simd
		for(int i = 1; i < N+1; i++)
		{
			for(int j = 1; j < N+1; j++)
			{
				FUNCTION_VECTOR(pos++) = Forcing(GRID(j),GRID(i));
			}
		}
		/*Initalize the boundaries*/
		matrix<T> LEFT(N,1),RIGHT(N,1),BOTTOM(N,1),TOP(N,1);
#pragma omp for simd
		for(int i = 0; i < N; i++)
		{
			LEFT(i) = boundary(0)(GRID(i+1));
			RIGHT(i) = boundary(1)(GRID(i+1));
			TOP(i) = boundary(2)(GRID(i+1));
			BOTTOM(i) = boundary(3)(GRID(i+1));
		}
		
		/*Generate Top part of the matrix (Bottom grid)*/
		unsigned int e = 0;
		for(int i = 0; i < N; i++)
		{
			unsigned int q0 = i;
			unsigned int q1 = i+N;
			if(i == 0)
			{
				if(boundary_type(0) == 'n')
				{
					if(boundary_type(3) == 'n')
					{
						MATRIX(i,q0) = (T)-2.;
						MATRIX(i,q0+1) = (T)1.;
						MATRIX(i,q1) = (T)1.;
						BOUNDARY_CORRECTION(i) = h*(LEFT(0) + BOTTOM(e));
					}
					else if(boundary_type(3) == 'd')
					{
						MATRIX(i,q0) = (T)-3.;
						MATRIX(i,q0+1) = (T)1.;
						MATRIX(i,q1) = (T)1.;
						BOUNDARY_CORRECTION(i) = h*BOTTOM(e) + -LEFT(0);
					}
					else {std::cout << "Improper type of conditions\nExiting\n"; exit(1);}
				}
				else if(boundary_type(0) == 'd')
				{
					if(boundary_type(3) == 'n')
					{
						MATRIX(i,q0) = (T)-3.;
						MATRIX(i,q0+1) = (T)1.;
						MATRIX(i,q1) = (T)1.;
						BOUNDARY_CORRECTION(i) = h*BOTTOM(e) + -LEFT(0);
					}
					else if(boundary_type(3) == 'd')
					{
						MATRIX(i,q0) = (T)-4.;
						MATRIX(i,q0+1) = (T)1.;
						MATRIX(i,q1) = (T)1.;
						BOUNDARY_CORRECTION(i) = -LEFT(0) + -BOTTOM(e);
					}
				}
				else {std::cout << "Improper type of conditions\nExiting\n"; exit(1);}
			}
			else if(i == N-1)
			{
				if(boundary_type(1) == 'n')
				{
					if(boundary_type(3) == 'n')
					{
						MATRIX(i,q0) = (T)-2.;
						MATRIX(i,q0-1) = (T)1.;
						MATRIX(i,q1) = (T)1.;
						BOUNDARY_CORRECTION(i) = h*(-RIGHT(0) + BOTTOM(e));
					}
					else if(boundary_type(3) == 'd')
					{
						MATRIX(i,q0) = (T)-3.;
						MATRIX(i,q0-1) = (T)1.;
						MATRIX(i,q1) = (T)1.;
						BOUNDARY_CORRECTION(i) = h*BOTTOM(e) + -RIGHT(0);
					}
					else {std::cout << "Improper type of conditions\nExiting\n"; exit(1);}
				}
				else if(boundary_type(1) == 'd')
				{
					if(boundary_type(3) == 'n')
					{
						MATRIX(i,q0) = (T)-3.;
						MATRIX(i,q0-1) = (T)1.;
						MATRIX(i,q1) = (T)1.;
						BOUNDARY_CORRECTION(i) = h*BOTTOM(e) + -RIGHT(0);
					}
					else if(boundary_type(3) == 'd')
					{
						MATRIX(i,q0) = (T)-4.;
						MATRIX(i,q0-1) = (T)1.;
						MATRIX(i,q1) = (T)1.;
						BOUNDARY_CORRECTION(i) = -RIGHT(0) + -BOTTOM(e);
					}
				}
				else {std::cout << "Improper type of conditions\nExiting\n"; exit(1);}
			}
			else
			{
				if(boundary_type(3) == 'n')
				{
					MATRIX(i,q0) = (T)-3.;
					MATRIX(i,q0+1) = (T)1.;
					MATRIX(i,q0-1) = (T)1.;
					MATRIX(i,q1) = (T)1.;
					BOUNDARY_CORRECTION(i)= h*BOTTOM(e);
				}
				else if(boundary_type(3) == 'd')
				{
					MATRIX(i,q0) = (T)-4.;
					MATRIX(i,q0+1) = (T)1.;
					MATRIX(i,q0-1) = (T)1.;
					MATRIX(i,q1) = (T)1.;
					BOUNDARY_CORRECTION(i)= -BOTTOM(e);
				}
				else {std::cout << "Improper type of conditions\nExiting\n"; exit(1);}
			}
			e++;
		}
			
		/*Generate bottom part of the matrix(top of grid)*/
		e = 0;
		for(int i = N*N - N; i < N*N; i++)
		{
			unsigned int q0 = i;
			unsigned int q1 = i-N;
			if(i == N*N - N)
			{
				if(boundary_type(0) == 'n')
				{
					if(boundary_type(2) == 'n')
					{
						MATRIX(i,q0) = (T)-2.;
						MATRIX(i,q0+1) = (T)1.;
						MATRIX(i,q1) = (T)1.;
						BOUNDARY_CORRECTION(i) = h*(LEFT(N-1) + -TOP(e));
					}
					else if(boundary_type(2) == 'd')
					{
						MATRIX(i,q0) = (T)-3.;
						MATRIX(i,q0+1) = (T)1.;
						MATRIX(i,q1) = (T)1.;
						BOUNDARY_CORRECTION(i) = -h*TOP(e) + -LEFT(N-1);
					}
					else {std::cout << "Improper type of conditions\nExiting\n"; exit(1);}
				}
				else if(boundary_type(0) == 'd')
				{
					if(boundary_type(2) == 'n')
					{
						MATRIX(i,q0) = (T)-3.;
						MATRIX(i,q0+1) = (T)1.;
						MATRIX(i,q1) = (T)1.;
						BOUNDARY_CORRECTION(i) = -h*TOP(e) + -LEFT(N-1);
					}
					else if(boundary_type(2) == 'd')
					{
						MATRIX(i,q0) = (T)-4.;
						MATRIX(i,q0+1) = (T)1.;
						MATRIX(i,q1) = (T)1.;
						BOUNDARY_CORRECTION(i) = -LEFT(N-1) + -TOP(e);
					}
				}
				else {std::cout << "Improper type of conditions\nExiting\n"; exit(1);}
			}
			else if(i == N*N-1)
			{
				if(boundary_type(1) == 'n')
				{
					if(boundary_type(2) == 'n')
					{
						MATRIX(i,q0) = (T)-2.;
						MATRIX(i,q0-1) = (T)1.;
						MATRIX(i,q1) = (T)1.;
						BOUNDARY_CORRECTION(i) = -h*(RIGHT(N-1) + TOP(e));
					}
					else if(boundary_type(2) == 'd')
					{
						MATRIX(i,q0) = (T)-3.;
						MATRIX(i,q0-1) = (T)1.;
						MATRIX(i,q1) = (T)1.;
						BOUNDARY_CORRECTION(i) = -h*TOP(e) + -RIGHT(N-1);
					}
					else {std::cout << "Improper type of conditions\nExiting\n"; exit(1);}
				}
				else if(boundary_type(1) == 'd')
				{
					if(boundary_type(2) == 'n')
					{
						MATRIX(i,q0) = (T)-3.;
						MATRIX(i,q0-1) = (T)1.;
						MATRIX(i,q1) = (T)1.;
						BOUNDARY_CORRECTION(i) = -h*TOP(e) + -RIGHT(N-1);
					}
					else if(boundary_type(2) == 'd')
					{
						MATRIX(i,q0) = (T)-4.;
						MATRIX(i,q0-1) = (T)1.;
						MATRIX(i,q1) = (T)1.;
						BOUNDARY_CORRECTION(i) = -RIGHT(N-1) + -TOP(e);
					}
				}
				else {std::cout << "Improper type of conditions\nExiting\n"; exit(1);}
			}
			else
			{
				if(boundary_type(2) == 'n')
				{
					MATRIX(i,q0) = (T)-3.;
					MATRIX(i,q0+1) = (T)1.;
					MATRIX(i,q0-1) = (T)1.;
					MATRIX(i,q1) = (T)1.;
					BOUNDARY_CORRECTION(i)= -h*TOP(e);
				}
				else if(boundary_type(2) == 'd')
				{
					MATRIX(i,q0) = (T)-4.;
					MATRIX(i,q0+1) = (T)1.;
					MATRIX(i,q0-1) = (T)1.;
					MATRIX(i,q1) = (T)1.;
					BOUNDARY_CORRECTION(i)= -TOP(e);
				}
				else {std::cout << "Improper type of conditions\nExiting\n"; exit(1);}
			}
			e++;
		}
			
		/*Generate internal part of the matrix*/
		unsigned int r = 1;
		unsigned int s = 1;
		for(int i = 0; i < N; i++)
		{
			for(int j = 1; j < N-1; j++)
			{
				unsigned int q = i + j*N;
				if(i == 0)
				{
					if(boundary_type(0) == 'n')
					{
						MATRIX(q,q) = (T)-3.;
						MATRIX(q,q+1) = (T)1.;
						MATRIX(q,q+N) = (T)1.;
						MATRIX(q,q-N) = (T)1.;
						BOUNDARY_CORRECTION(q) = h*LEFT(r);
					}
					else if(boundary_type(0) == 'd')
					{
						MATRIX(q,q) = (T)-4.;
						MATRIX(q,q+1) = (T)1.;
						MATRIX(q,q+N) = (T)1.;
						MATRIX(q,q-N) = (T)1.;
						BOUNDARY_CORRECTION(q) = -LEFT(r);
					}
					else {std::cout << "Improper type of conditions\nExiting\n"; exit(1);}
					r++;
				}
				else if(i == N-1)
				{
					if(boundary_type(1) == 'n')
					{
						MATRIX(q,q) = (T)-3.;
						MATRIX(q,q-1) = (T)1.;
						MATRIX(q,q+N) = (T)1.;
						MATRIX(q,q-N) = (T)1.;
						BOUNDARY_CORRECTION(q) = -h*RIGHT(s);
					}
					else if(boundary_type(1) == 'd')
					{
						MATRIX(q,q) = (T)-4.;
						MATRIX(q,q-1) = (T)1.;
						MATRIX(q,q+N) = (T)1.;
						MATRIX(q,q-N) = (T)1.;
						BOUNDARY_CORRECTION(q) = -RIGHT(s);
					}
					else {std::cout << "Improper type of conditions\nExiting\n"; exit(1);}
					s++;
				}
				else
				{
					MATRIX(q,q) = (T)-4.;
					MATRIX(q,q-1) = (T)1.;
					MATRIX(q,q+1) = (T)1.;
					MATRIX(q,q+N) = (T)1.;
					MATRIX(q,q-N) = (T)1.;
				}
			}	
		}
		std::cout << "Matrix and Vectors Initalized\n";
		//std::cout << MATRIX << FUNCTION_VECTOR << BOUNDARY_CORRECTION;
	}
	
	/*LU SOLVER!*/
	template<class T>
	matrix<T> laplace<T>::lu_solve(void)
	{
		auto t1 = std::chrono::high_resolution_clock::now();
		const unsigned int N = GRID.numel()-2;
		matrix<T> sol = gauss<T>::solve(MATRIX,BOUNDARY_CORRECTION - h*h*FUNCTION_VECTOR,piviot_args::nopiviot); std::cout << "SOLVED!\n";
		//std::cout << sol;
		/*Save the results to a file*/
		FILE *fin = fopen("Laplace_sol.txt","w"); if(!fin){std::cout << "Runtime Error\nCannot open file\nExiting\n";exit(1);} else {}
		for(int i = 0; i < N; i++)
		{
			for(int j = 0; j < N; j++)
			{
				fprintf(fin,"%.14f\t%.14f\t%.14f\n",GRID(j+1),GRID(i+1),sol(i*N + j));
			}
		}
		/*Plot boundary conditions*/
		matrix<T> GRID_p = GRID.get_row_range(0,1,N);
		/*bottom*/
		for(int i = 0; i < N; i++)
		{
			if(boundary_type(3) == 'd')
			{
				fprintf(fin,"%.14f\t%.14f\t%.14f\n",GRID_p(i),GRID(0),boundary(3)(GRID_p(i)));
			}
			else if(boundary_type(3) == 'n')
			{
				fprintf(fin,"%.14f\t%.14f\t%.14f\n",GRID_p(i),GRID(0),sol(i) - h*boundary(3)(GRID_p(i)));
			}
		}
		/*Top*/
		for(int i = 0; i < N; i++)
		{
			if(boundary_type(2) == 'd')
			{
				fprintf(fin,"%.14f\t%.14f\t%.14f\n",GRID_p(i),GRID(N+1),boundary(2)(GRID_p(i)));
			}
			else if(boundary_type(2) == 'n')
			{
				fprintf(fin,"%.14f\t%.14f\t%.14f\n",GRID_p(i),GRID(N+1),sol(N*N - N + i) + h*boundary(2)(GRID_p(i)));
			}
		}
		/*Left*/
		for(int i = 0; i < N; i++)
		{
			if(boundary_type(0) == 'd')
			{
				fprintf(fin,"%.14f\t%.14f\t%.14f\n",GRID(0),GRID_p(i),boundary(0)(GRID_p(i)));
			}
			else if(boundary_type(0) == 'n')
			{
				fprintf(fin,"%.14f\t%.14f\t%.14f\n",GRID(0),GRID_p(i),sol(i*N) - h*boundary(0)(GRID_p(i)));
			}
		}
		/*Right*/
		for(int i = 0; i < N; i++)
		{
			if(boundary_type(1) == 'd')
			{
				fprintf(fin,"%.14f\t%.14f\t%.14f\n",GRID(N+1),GRID_p(i),boundary(1)(GRID_p(i)));
			}
			else if(boundary_type(1) == 'n')
			{
				fprintf(fin,"%.14f\t%.14f\t%.14f\n",GRID(N+1),GRID_p(i),sol(i*N + N - 1) + h*boundary(1)(GRID_p(i)));
			}
		}
		fclose(fin);
		std::cout << "Solved in " << std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - t1).count() << " seconds\n";
		return sol;
	}
	
	template<class T>
	matrix<T> laplace<T>::j_solve(unsigned int max_iter,T tolerence)
	{
		auto t1 = std::chrono::high_resolution_clock::now();
		const unsigned int N = GRID.numel()-2;
		matrix<T> sol = jacobi<T>::solve(MATRIX,BOUNDARY_CORRECTION - h*h*FUNCTION_VECTOR,tolerence,max_iter); std::cout << "SOLVED!\n";
		/*Save the results to a file*/
		FILE *fin = fopen("Laplace_sol.txt","w"); if(!fin){std::cout << "Runtime Error\nCannot open file\nExiting\n";exit(1);} else {}
		for(int i = 0; i < N; i++)
		{
			for(int j = 0; j < N; j++)
			{
				fprintf(fin,"%.14f\t%.14f\t%.14f\n",GRID(j+1),GRID(i+1),sol(i*N + j));
			}
		}
		
		/*Plot boundary conditions*/
		matrix<T> GRID_p = GRID.get_row_range(0,1,N);
		/*bottom*/
		for(int i = 0; i < N; i++)
		{
			if(boundary_type(3) == 'd')
			{
				fprintf(fin,"%.14f\t%.14f\t%.14f\n",GRID_p(i),GRID(0),boundary(3)(GRID_p(i)));
			}
			else if(boundary_type(3) == 'n')
			{
				fprintf(fin,"%.14f\t%.14f\t%.14f\n",GRID_p(i),GRID(0),sol(i) - h*boundary(3)(GRID_p(i)));
			}
		}
		/*Top*/
		for(int i = 0; i < N; i++)
		{
			if(boundary_type(2) == 'd')
			{
				fprintf(fin,"%.14f\t%.14f\t%.14f\n",GRID_p(i),GRID(N+1),boundary(2)(GRID_p(i)));
			}
			else if(boundary_type(2) == 'n')
			{
				fprintf(fin,"%.14f\t%.14f\t%.14f\n",GRID_p(i),GRID(N+1),sol(N*N - N + i) + h*boundary(2)(GRID_p(i)));
			}
		}
		/*Left*/
		for(int i = 0; i < N; i++)
		{
			if(boundary_type(0) == 'd')
			{
				fprintf(fin,"%.14f\t%.14f\t%.14f\n",GRID(0),GRID_p(i),boundary(0)(GRID_p(i)));
			}
			else if(boundary_type(0) == 'n')
			{
				fprintf(fin,"%.14f\t%.14f\t%.14f\n",GRID(0),GRID_p(i),sol(i*N) - h*boundary(0)(GRID_p(i)));
			}
		}
		/*Right*/
		for(int i = 0; i < N; i++)
		{
			if(boundary_type(1) == 'd')
			{
				fprintf(fin,"%.14f\t%.14f\t%.14f\n",GRID(N+1),GRID_p(i),boundary(1)(GRID_p(i)));
			}
			else if(boundary_type(1) == 'n')
			{
				fprintf(fin,"%.14f\t%.14f\t%.14f\n",GRID(N+1),GRID_p(i),sol(i*N + N - 1) + h*boundary(1)(GRID_p(i)));
			}
		}
		
		fclose(fin);
		std::cout << "Solved in " << std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - t1).count() << " seconds\n";
		return sol;
	}
	
	/*//////////////////////////////////////////////////////////////////////////////////////////////*/
	/*For High Precision*/
	template<>
	matrix<long double> laplace<long double>::lu_solve(void)
	{
		auto t1 = std::chrono::high_resolution_clock::now();
		const unsigned int N = GRID.numel()-2;
		matrix<long double> sol = gauss<long double>::solve(MATRIX,BOUNDARY_CORRECTION - h*h*FUNCTION_VECTOR,piviot_args::nopiviot); std::cout << "SOLVED!\n";
		//std::cout << sol;
		/*Save the results to a file*/
		FILE *fin = fopen("Laplace_sol.txt","w"); if(!fin){std::cout << "Runtime Error\nCannot open file\nExiting\n";exit(1);} else {}
		for(int i = 0; i < N; i++)
		{
			for(int j = 0; j < N; j++)
			{
				fprintf(fin,"%.14Lf\t%.14Lf\t%.14Lf\n",GRID(j+1),GRID(i+1),sol(i*N + j));
			}
		}
		/*Plot boundary conditions*/
		matrix<long double> GRID_p = GRID.get_row_range(0,1,N);
		/*bottom*/
		for(int i = 0; i < N; i++)
		{
			if(boundary_type(3) == 'd')
			{
				fprintf(fin,"%.14Lf\t%.14Lf\t%.14Lf\n",GRID_p(i),GRID(0),boundary(3)(GRID_p(i)));
			}
			else if(boundary_type(3) == 'n')
			{
				fprintf(fin,"%.14Lf\t%.14Lf\t%.14Lf\n",GRID_p(i),GRID(0),sol(i) - h*boundary(3)(GRID_p(i)));
			}
		}
		/*Top*/
		for(int i = 0; i < N; i++)
		{
			if(boundary_type(2) == 'd')
			{
				fprintf(fin,"%.14Lf\t%.14Lf\t%.14Lf\n",GRID_p(i),GRID(N+1),boundary(2)(GRID_p(i)));
			}
			else if(boundary_type(2) == 'n')
			{
				fprintf(fin,"%.14Lf\t%.14Lf\t%.14Lf\n",GRID_p(i),GRID(N+1),sol(N*N - N + i) + h*boundary(2)(GRID_p(i)));
			}
		}
		/*Left*/
		for(int i = 0; i < N; i++)
		{
			if(boundary_type(0) == 'd')
			{
				fprintf(fin,"%.14Lf\t%.14Lf\t%.14Lf\n",GRID(0),GRID_p(i),boundary(0)(GRID_p(i)));
			}
			else if(boundary_type(0) == 'n')
			{
				fprintf(fin,"%.14Lf\t%.14Lf\t%.14Lf\n",GRID(0),GRID_p(i),sol(i*N) - h*boundary(0)(GRID_p(i)));
			}
		}
		/*Right*/
		for(int i = 0; i < N; i++)
		{
			if(boundary_type(1) == 'd')
			{
				fprintf(fin,"%.14Lf\t%.14Lf\t%.14Lf\n",GRID(N+1),GRID_p(i),boundary(1)(GRID_p(i)));
			}
			else if(boundary_type(1) == 'n')
			{
				fprintf(fin,"%.14Lf\t%.14Lf\t%.14Lf\n",GRID(N+1),GRID_p(i),sol(i*N + N - 1) + h*boundary(1)(GRID_p(i)));
			}
		}
		fclose(fin);
		std::cout << "Solved in " << std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - t1).count() << " seconds\n";
		return sol;
	}
	
	template<>
	matrix<long double> laplace<long double>::j_solve(unsigned int max_iter,long double tolerence)
	{
		auto t1 = std::chrono::high_resolution_clock::now();
		const unsigned int N = GRID.numel()-2;
		matrix<long double> sol = jacobi<long double>::solve(MATRIX,BOUNDARY_CORRECTION - h*h*FUNCTION_VECTOR,tolerence,max_iter); std::cout << "SOLVED!\n";
		/*Save the results to a file*/
		FILE *fin = fopen("Laplace_sol.txt","w"); if(!fin){std::cout << "Runtime Error\nCannot open file\nExiting\n";exit(1);} else {}
		for(int i = 0; i < N; i++)
		{
			for(int j = 0; j < N; j++)
			{
				fprintf(fin,"%.14Lf\t%.14Lf\t%.14Lf\n",GRID(j+1),GRID(i+1),sol(i*N + j));
			}
		}
		
		/*Plot boundary conditions*/
		matrix<long double> GRID_p = GRID.get_row_range(0,1,N);
		/*bottom*/
		for(int i = 0; i < N; i++)
		{
			if(boundary_type(3) == 'd')
			{
				fprintf(fin,"%.14Lf\t%.14Lf\t%.14Lf\n",GRID_p(i),GRID(0),boundary(3)(GRID_p(i)));
			}
			else if(boundary_type(3) == 'n')
			{
				fprintf(fin,"%.14Lf\t%.14Lf\t%.14Lf\n",GRID_p(i),GRID(0),sol(i) - h*boundary(3)(GRID_p(i)));
			}
		}
		/*Top*/
		for(int i = 0; i < N; i++)
		{
			if(boundary_type(2) == 'd')
			{
				fprintf(fin,"%.14Lf\t%.14Lf\t%.14Lf\n",GRID_p(i),GRID(N+1),boundary(2)(GRID_p(i)));
			}
			else if(boundary_type(2) == 'n')
			{
				fprintf(fin,"%.14Lf\t%.14Lf\t%.14Lf\n",GRID_p(i),GRID(N+1),sol(N*N - N + i) + h*boundary(2)(GRID_p(i)));
			}
		}
		/*Left*/
		for(int i = 0; i < N; i++)
		{
			if(boundary_type(0) == 'd')
			{
				fprintf(fin,"%.14Lf\t%.14Lf\t%.14Lf\n",GRID(0),GRID_p(i),boundary(0)(GRID_p(i)));
			}
			else if(boundary_type(0) == 'n')
			{
				fprintf(fin,"%.14Lf\t%.14Lf\t%.14Lf\n",GRID(0),GRID_p(i),sol(i*N) - h*boundary(0)(GRID_p(i)));
			}
		}
		/*Right*/
		for(int i = 0; i < N; i++)
		{
			if(boundary_type(1) == 'd')
			{
				fprintf(fin,"%.14Lf\t%.14Lf\t%.14Lf\n",GRID(N+1),GRID_p(i),boundary(1)(GRID_p(i)));
			}
			else if(boundary_type(1) == 'n')
			{
				fprintf(fin,"%.14Lf\t%.14Lf\t%.14Lf\n",GRID(N+1),GRID_p(i),sol(i*N + N - 1) + h*boundary(1)(GRID_p(i)));
			}
		}
		
		fclose(fin);
		std::cout << "Solved in " << std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - t1).count() << " seconds\n";
		return sol;
	}
}
