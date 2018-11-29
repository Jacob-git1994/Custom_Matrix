#pragma once
#include "matrix.h"
#include <iostream>
#include <stdlib.h>

namespace nicole
{
	template<class T>
	class conjgrad
	{
	public:		
		static matrix<T> solve(const matrix<T>&,const matrix<T>&);
		static matrix<T> solve(const matrix<T>&,const matrix<T>&,T);
		static matrix<T> solve(const matrix<T>&,const matrix<T>&,unsigned int);
		static matrix<T> solve(const matrix<T>&,const matrix<T>&,const T,const unsigned int);
		
		static matrix<T> solve(const matrix<T>&,const matrix<T>&,matrix<T>&);
		static matrix<T> solve(const matrix<T>&,const matrix<T>&,matrix<T>&,const T);
		static matrix<T> solve(const matrix<T>&,const matrix<T>&,matrix<T>&,const unsigned int);
		static matrix<T> solve(const matrix<T>&,const matrix<T>&,matrix<T>&,const T,const unsigned int);
	};

	template<class T>
	matrix<T> conjgrad<T>::solve(const matrix<T> &A,const matrix<T> &b)
	{
		/*Check to see if A is square and no column mismatch*/
		if(!A.is_square() || A.cols() != b.rows() || b.cols() != 1) {std::cout << "Runtime Error\nMatrix not the same size or not square\n"; exit(1);} 
		else {}
		
		/*Initalize the tolerences and max iterations*/
		const T tol = .000001;
		const unsigned int MAX_ITER = b.numel();
		unsigned int count = 0;
		
		/*Initalize our alpha and beta*/
		T alpha = (T)0.;
		T beta = (T)0.;
		
		/*Initalize the residual and path direction*/
		matrix<T> x = matrix<T>::zeros(b.numel(),1);
		matrix<T> res = b - A*x;
		matrix<T> res_2 = res;
		matrix<T> p = res;
		
		/*Start iterations*/
		while((norm_2(res) >= tol) && (count < MAX_ITER))
		{
			alpha = (res.t()*res)(0)/(p.t()*A*p)(0);
			x += alpha*p;
			res_2 = res - alpha*A*p;
			/*Update search direction*/
			beta = (res_2.t()*res_2)(0)/(res.t()*res)(0);
			p = res_2 + beta*p;
			/*update results*/
			res = res_2;
			count++;
		}
		return x;
	}

	template<class T>
	matrix<T> conjgrad<T>::solve(const matrix<T> &A,const matrix<T> &b,const T tol)
	{
		/*Check to see if A is square and no column mismatch*/
		if(!A.is_square() || A.cols() != b.rows() || b.cols() != 1) {std::cout << "Runtime Error\nMatrix not the same size or not square\n"; exit(1);} 
		else {}
		
		/*Initalize the tolerences and max iterations*/
		const unsigned int MAX_ITER = b.numel();
		unsigned int count = 0;
		
		/*Initalize our alpha and beta*/
		T alpha = (T)0.;
		T beta = (T)0.;
		
		/*Initalize the residual and path direction*/
		matrix<T> x = matrix<T>::zeros(b.numel(),1);
		matrix<T> res = b - A*x;
		matrix<T> res_2 = res;
		matrix<T> p = res;
		
		/*Start iterations*/
		while((norm_2(res) >= tol) && (count < MAX_ITER))
		{
			alpha = (res.t()*res)(0)/(p.t()*A*p)(0);
			x += alpha*p;
			res_2 = res - alpha*A*p;
			/*Update search direction*/
			beta = (res_2.t()*res_2)(0)/(res.t()*res)(0);
			p = res_2 + beta*p;
			/*update results*/
			res = res_2;
			count++;
		}
		return x;
	}

	template<class T>
	matrix<T> conjgrad<T>::solve(const matrix<T> &A,const matrix<T> &b,const unsigned int MAX_ITER)
	{
		/*Check to see if A is square and no column mismatch*/
		if(!A.is_square() || A.cols() != b.rows() || b.cols() != 1) {std::cout << "Runtime Error\nMatrix not the same size or not square\n"; exit(1);} 
		else {}
		
		/*Initalize the tolerences and max iterations*/
		const T tol = .000001;
		unsigned int count = 0;
		
		/*Initalize our alpha and beta*/
		T alpha = (T)0.;
		T beta = (T)0.;
		
		/*Initalize the residual and path direction*/
		matrix<T> x = matrix<T>::zeros(b.numel(),1);
		matrix<T> res = b - A*x;
		matrix<T> res_2 = res;
		matrix<T> p = res;
		
		/*Start iterations*/
		while((norm_2(res) >= tol) && (count < MAX_ITER))
		{
			alpha = (res.t()*res)(0)/(p.t()*A*p)(0);
			x += alpha*p;
			res_2 = res - alpha*A*p;
			/*Update search direction*/
			beta = (res_2.t()*res_2)(0)/(res.t()*res)(0);
			p = res_2 + beta*p;
			/*update results*/
			res = res_2;
			count++;
		}
		return x;
	}

	template<class T>
	matrix<T> conjgrad<T>::solve(const matrix<T> &A,const matrix<T> &b,const T tol,const unsigned int MAX_ITER)
	{
		/*Check to see if A is square and no column mismatch*/
		if(!A.is_square() || A.cols() != b.rows() || b.cols() != 1) {std::cout << "Runtime Error\nMatrix not the same size or not square\n"; exit(1);} 
		else {}
		
		/*Initalize the tolerences and max iterations*/
		unsigned int count = 0;
		
		/*Initalize our alpha and beta*/
		T alpha = (T)0.;
		T beta = (T)0.;
		
		/*Initalize the residual and path direction*/
		matrix<T> x = matrix<T>::zeros(b.numel(),1);
		matrix<T> res = b - A*x;
		matrix<T> res_2 = res;
		matrix<T> p = res;
		
		/*Start iterations*/
		while((norm_2(res) >= tol) && (count < MAX_ITER))
		{
			alpha = (res.t()*res)(0)/(p.t()*A*p)(0);
			x += alpha*p;
			res_2 = res - alpha*A*p;
			/*Update search direction*/
			beta = (res_2.t()*res_2)(0)/(res.t()*res)(0);
			p = res_2 + beta*p;
			/*update results*/
			res = res_2;
			count++;
		}
		return x;
	}
	
	template<class T>
	matrix<T> conjgrad<T>::solve(const matrix<T> &A,const matrix<T> &b,matrix<T> &res_norm)
	{
		/*Check to see if A is square and no column mismatch*/
		if(!A.is_square() || A.cols() != b.rows() || b.cols() != 1) {std::cout << "Runtime Error\nMatrix not the same size or not square\n"; exit(1);} 
		else {}
		
		/*Initalize the tolerences and max iterations*/
		const T tol = .000001;
		const unsigned int MAX_ITER = b.numel();
		unsigned int count = 0;
		
		/*Initalize our alpha and beta*/
		T alpha = (T)0.;
		T beta = (T)0.;
		
		/*Initalize the residual and path direction*/
		matrix<T> x = matrix<T>::zeros(b.numel(),1);
		matrix<T> res = b - A*x;
		matrix<T> res_2 = res;
		matrix<T> p = res;
		
		/*Vector to save the residue in*/
		res_norm = matrix<T>::zeros(MAX_ITER,1);
			
		/*Start iterations*/
		while((norm_2(res) >= tol) && (count < MAX_ITER))
		{
			res_norm(count) = norm_2(res);
			alpha = (res.t()*res)(0)/(p.t()*A*p)(0);
			x += alpha*p;
			res_2 = res - alpha*A*p;
			/*Update search direction*/
			beta = (res_2.t()*res_2)(0)/(res.t()*res)(0);
			p = res_2 + beta*p;
			/*update results*/
			res = res_2;
			count++;
		}
		return x;
	}
	
	template<class T>
	matrix<T> conjgrad<T>::solve(const matrix<T> &A,const matrix<T> &b,matrix<T> &res_norm,const T tol)
	{
		/*Check to see if A is square and no column mismatch*/
		if(!A.is_square() || A.cols() != b.rows() || b.cols() != 1) {std::cout << "Runtime Error\nMatrix not the same size or not square\n"; exit(1);} 
		else {}
		
		/*Initalize the tolerences and max iterations*/
		const unsigned int MAX_ITER = b.numel();
		unsigned int count = 0;
		
		/*Initalize our alpha and beta*/
		T alpha = (T)0.;
		T beta = (T)0.;
		
		/*Initalize the residual and path direction*/
		matrix<T> x = matrix<T>::zeros(b.numel(),1);
		matrix<T> res = b - A*x;
		matrix<T> res_2 = res;
		matrix<T> p = res;
		
		/*Vector to save the residue in*/
		res_norm = matrix<T>::zeros(MAX_ITER,1);
			
		/*Start iterations*/
		while((norm_2(res) >= tol) && (count < MAX_ITER))
		{
			res_norm(count) = norm_2(res);
			alpha = (res.t()*res)(0)/(p.t()*A*p)(0);
			x += alpha*p;
			res_2 = res - alpha*A*p;
			/*Update search direction*/
			beta = (res_2.t()*res_2)(0)/(res.t()*res)(0);
			p = res_2 + beta*p;
			/*update results*/
			res = res_2;
			count++;
		}
		return x;
	}
	
	template<class T>
	matrix<T> conjgrad<T>::solve(const matrix<T> &A,const matrix<T> &b,matrix<T> &res_norm,const unsigned int MAX_ITER)
	{
		/*Check to see if A is square and no column mismatch*/
		if(!A.is_square() || A.cols() != b.rows() || b.cols() != 1) {std::cout << "Runtime Error\nMatrix not the same size or not square\n"; exit(1);} 
		else {}
		
		/*Initalize the tolerences and max iterations*/
		const T tol = .000001;
		unsigned int count = 0;
		
		/*Initalize our alpha and beta*/
		T alpha = (T)0.;
		T beta = (T)0.;
		
		/*Initalize the residual and path direction*/
		matrix<T> x = matrix<T>::zeros(b.numel(),1);
		matrix<T> res = b - A*x;
		matrix<T> res_2 = res;
		matrix<T> p = res;
		
		/*Vector to save the residue in*/
		res_norm = matrix<T>::zeros(MAX_ITER,1);
			
		/*Start iterations*/
		while((norm_2(res) >= tol) && (count < MAX_ITER))
		{
			res_norm(count) = norm_2(res);
			alpha = (res.t()*res)(0)/(p.t()*A*p)(0);
			x += alpha*p;
			res_2 = res - alpha*A*p;
			/*Update search direction*/
			beta = (res_2.t()*res_2)(0)/(res.t()*res)(0);
			p = res_2 + beta*p;
			/*update results*/
			res = res_2;
			count++;
		}
		return x;
	}
	
	template<class T>
	matrix<T> conjgrad<T>::solve(const matrix<T> &A,const matrix<T> &b,matrix<T> &res_norm,const T tol,const unsigned int MAX_ITER)
	{
		/*Check to see if A is square and no column mismatch*/
		if(!A.is_square() || A.cols() != b.rows() || b.cols() != 1) {std::cout << "Runtime Error\nMatrix not the same size or not square\n"; exit(1);} 
		else {}
		
		/*Initalize the tolerences and max iterations*/
		unsigned int count = 0;
		
		/*Initalize our alpha and beta*/
		T alpha = (T)0.;
		T beta = (T)0.;
		
		/*Initalize the residual and path direction*/
		matrix<T> x = matrix<T>::zeros(b.numel(),1);
		matrix<T> res = b - A*x;
		matrix<T> res_2 = res;
		matrix<T> p = res;
		
		/*Vector to save the residue in*/
		res_norm = matrix<T>::zeros(MAX_ITER,1);
			
		/*Start iterations*/
		while((norm_2(res) >= tol) && (count < MAX_ITER))
		{
			res_norm(count) = norm_2(res);
			alpha = (res.t()*res)(0)/(p.t()*A*p)(0);
			x += alpha*p;
			res_2 = res - alpha*A*p;
			/*Update search direction*/
			beta = (res_2.t()*res_2)(0)/(res.t()*res)(0);
			p = res_2 + beta*p;
			/*update results*/
			res = res_2;
			count++;
		}
		return x;
	}
}