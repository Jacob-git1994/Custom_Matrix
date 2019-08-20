#pragma once
#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <iomanip>
#include <complex>

namespace nicole
{	
	static const long double pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286;
	static const long double e = std::exp(1);
	
	template<class T>
	class matrix
	{
	private:
		unsigned int row,col;
		unsigned int if_in(unsigned int);
		unsigned int if_in(unsigned int,unsigned int);
		T *data;
	public:
		/*Constructors*/
		matrix(void);
		matrix(unsigned int); /*Square matrix*/
		matrix(unsigned int,unsigned int); /*Different sized matrix*/
		matrix(unsigned int,const T[]); /*Copies elements in col order*/
		matrix(unsigned int,unsigned int, const T[]); /*Copies elements in col order for different sizes*/
		matrix(const matrix<T>&); /*Copy constructor*/
		matrix(const T&);
		/*Destructor*/
		~matrix(void);
		/*Printing out data*/
		template<class U> friend std::ostream & operator<<(std::ostream&,const matrix<U>&);
		/*Indexing*/
		T & operator()(unsigned int);
		T operator()(unsigned int) const;
		T operator()(unsigned int,unsigned int) const;
		T & operator()(unsigned int,unsigned int);
		/*Read only indexing ->used for const matrix references*/
		T at(unsigned int) const;
		T at(unsigned int,unsigned int) const;
		/*Arethmetic operations*/
		matrix<T> & operator=(const matrix<T>&);
		matrix<T> & operator=(const T&);
		matrix<T> & operator+(void);
		matrix<T> & operator+=(const matrix<T>&);
		matrix<T> & operator-=(const matrix<T>&);
		matrix<T> & operator*=(T);
		matrix<T> & operator/=(T);
		/*Element operations*/
		matrix<T> & operator&&(const matrix<T>&);
		matrix<T> & operator||(const matrix<T>&);
		/*Friend Arethmetic operations*/
		template<class U> friend matrix<U> operator-(const matrix<U>&);
		template<class U> friend matrix<U> operator-(const matrix<U>&,const matrix<U>&);
		template<class U> friend matrix<U> operator+(const matrix<U>&,const matrix<U>&);
		template<class U> friend matrix<U> operator*(U,const matrix<U>&);
		template<class U> friend matrix<U> operator*(const matrix<U>&,U);
		template<class U> friend matrix<U> operator/(U,const matrix<U>&);
		template<class U> friend matrix<U> operator/(const matrix<U>&,U);
		/*Friend Element operations*/
		template<class U> friend matrix<U> operator&(const matrix<U>&,const matrix<U>&);
		template<class U> friend matrix<U> operator|(const matrix<U>&,const matrix<U>&);
		/*Generating friend functions*/
		static matrix<T> zeros(unsigned int);
		static matrix<T> zeros(unsigned int,unsigned int);
		static matrix<T> ones(unsigned int);
		static matrix<T> ones(unsigned int,unsigned int);
		static matrix<T> randu(unsigned int);
		static matrix<T> randu(unsigned int,unsigned int);
		static matrix<T> randu(double,double,unsigned int);
		static matrix<T> randu(double,double,unsigned int,unsigned int);
		static matrix<T> eye(unsigned int);
		static matrix<T> linspace(T,T,unsigned int);
		static matrix<T> delspace(T,T,T);
		/*Control the random generation*/
		static void set_rand(void);
		static void set_rand(long long);
		/*Uninitalized matrix*/
		static matrix<T> shape(unsigned int);
		static matrix<T> shape(unsigned int,unsigned int);
		/*Normal random variables*/
		static matrix<T> randn(unsigned int);
		static matrix<T> randn(unsigned int,unsigned int);
		static matrix<T> randn(double,double,unsigned int);
		static matrix<T> randn(double,double,unsigned int,unsigned int);
		/*Matrix - Vector Operations*/
		matrix<T> t(void); /*Transpose*/
		matrix<T> t(void) const; /*Constant transpose*/
		template<class U> friend matrix<U> row_concat(const matrix<U>&,const matrix<U>&);
		template<class U> friend matrix<U> col_concat(const matrix<U>&,const matrix<U>&);
		unsigned int rows(void) const;
		unsigned int cols(void) const;
		unsigned int numel(void) const;
		matrix<T> vectorize(void);
		matrix<T> vectorize(unsigned int);
		matrix<T> & fill(T);
		/*Operations of chunks of vector or matrix >> doesnt support autominius functions*/
		matrix<T> & for_col(unsigned int,T (*)(T));
		matrix<T> & for_row(unsigned int,T (*)(T));
		matrix<T> & for_col_range(unsigned int,unsigned int,unsigned int,T (*)(T));
		matrix<T> & for_row_range(unsigned int,unsigned int,unsigned int,T (*)(T));
		matrix<T> & for_all(T (*)(T));
		/*Special method that returns the raw data*/
		T * raw_data(void) const {return data;}
		/*Swap rows or cols*/
		matrix<T> & row_swap(unsigned int,unsigned int);
		matrix<T> & col_swap(unsigned int,unsigned int);
		/*Returns the row,col as a copy as well as getting the diagonal (overload a friend operator for generating a diagonal)*/
		//matrix<T> get_row(unsigned int);
		matrix<T> get_row(unsigned int);
		matrix<T> get_col(unsigned int);
		matrix<T> get_row_range(unsigned int,unsigned int,unsigned int);
		matrix<T> get_col_range(unsigned int,unsigned int,unsigned int);
		matrix<T> diag(void);
		unsigned int is_in(T);
		unsigned int is_square(void) const;
		unsigned int is_vector(void) const;
		unsigned int is_row_vector(void) const;
		unsigned int is_col_vector(void) const;
		/*matrix multiplication*/
		template<class U> friend matrix<U> operator*(const matrix<U>&,const matrix<U>&);
		/*Banded matrix - Uses band on the left hand side*/
		template<class U> friend matrix<U> mat_mull_vec_band(const matrix<U>&,const matrix<U>&,unsigned int);
		/*Stat functions - more to come*/
		template<class U> friend U sum(const matrix<U>&);
		template<class U> friend U mean(const matrix<U>&);
		template<class U> friend U var(const matrix<U>&);
		template<class U> friend U sd(const matrix<U>&);
		template<class U> friend matrix<U> cov(const matrix<U>&,const matrix<U>&);
		template<class U> friend matrix<U> cov(const matrix<U>&);
		template<class U> friend U max(const matrix<U>&);
		template<class U> friend U min(const matrix<U>&);
		/*Functions to generate permutations forming matrix obects*/
		static matrix<T> randperm(unsigned int);
		static matrix<T> randperm(unsigned int,unsigned int);
 		/*General functions*/
		template<class U> friend U inf_norm(const matrix<U>&);
		template<class U> friend U one_norm(const matrix<U>&);
		template<class U> friend U norm_2(const matrix<U>&);
		template<class U> friend U p_norm(const matrix<U>&,double);
		template<class U> friend U rcond(const matrix<U>&);
		/*cmath functions overloads*/
		template<class U> friend matrix<U> cos(const matrix<U>&);
		template<class U> friend matrix<U> sin(const matrix<U>&);
		template<class U> friend matrix<U> tan(const matrix<U>&);
		template<class U> friend matrix<U> acos(const matrix<U>&);
		template<class U> friend matrix<U> asin(const matrix<U>&);
		template<class U> friend matrix<U> atan(const matrix<U>&);
		template<class U> friend matrix<U> cosh(const matrix<U>&);
		template<class U> friend matrix<U> sinh(const matrix<U>&);
		template<class U> friend matrix<U> tanh(const matrix<U>&);
		template<class U> friend matrix<U> acosh(const matrix<U>&);
		template<class U> friend matrix<U> asinh(const matrix<U>&);
		template<class U> friend matrix<U> atanh(const matrix<U>&);
		template<class U> friend matrix<U> exp(const matrix<U>&);
		template<class U> friend matrix<U> log(const matrix<U>&);
		template<class U> friend matrix<U> log10(const matrix<U>&);
		template<class U> friend matrix<U> log2(const matrix<U>&);
		template<class U> friend matrix<U> operator^(const matrix<U>&,double);
		template<class U> friend matrix<U> operator^(double,const matrix<U>&);
		template<class U> friend matrix<U> ceil(const matrix<U>&);
		template<class U> friend matrix<U> floor(const matrix<U>&);
		template<class U> friend matrix<U> abs(const matrix<U>&);
		template<class U> friend matrix<U> diff(const matrix<U>&); /*Calculates the diff in row wise format*/
		/*operator overloading to match matlab => allows quick operations to whole matrix*/
		template<class U> friend matrix<U> operator+(const matrix<U>&,U);
		template<class U> friend matrix<U> operator+(U,const matrix<U>&);
		template<class U> friend matrix<U> operator-(const matrix<U>&,U);
		template<class U> friend matrix<U> operator-(U,const matrix<U>&);
		matrix<T> & operator++(void);
		/*Logical operators*/
		template<class U> friend matrix<unsigned int> operator>(const matrix<U>&,U);
		template<class U> friend matrix<unsigned int> operator>(U,const matrix<U>&);
		template<class U> friend matrix<unsigned int> operator<(const matrix<U>&,U);
		template<class U> friend matrix<unsigned int> operator<(U,const matrix<U>&);
		template<class U> friend matrix<unsigned int> operator>=(const matrix<U>&,U);
		template<class U> friend matrix<unsigned int> operator>=(U,const matrix<U>&);
		template<class U> friend matrix<unsigned int> operator<=(const matrix<U>&,U);
		template<class U> friend matrix<unsigned int> operator<=(U,const matrix<U>&);
		template<class U> friend matrix<unsigned int> operator==(const matrix<U>&,U);
		template<class U> friend matrix<unsigned int> operator==(U,const matrix<U>&);
		template<class U> friend matrix<unsigned int> operator!=(const matrix<U>&,U);
		template<class U> friend matrix<unsigned int> operator!=(U,const matrix<U>&);
		/*Overload index operator to refind matrix with logical conditions*/
		matrix<T> operator()(const matrix<unsigned int>&); /*Vectorizes due to some issues when reconstructing matrix*/
		/*Turns the matrix into a matrix with index values*/
		matrix<unsigned int> find(const matrix<unsigned int>&);
		/*Set rows/cols/chuncks of the matrix*/
		matrix<T> & set_col(unsigned int,const matrix<T>&);
		matrix<T> & set_row(unsigned int,const matrix<T>&);
		
		/*ADD IN THIS CODE if needed*/
		///////////////////////////////////////////////////////////////////////////////////////////////////////////
		/*Set rows/cols/chuncks of the matrix*/
		matrix<T> & set_col_range(unsigned int,unsigned int,unsigned int,const matrix<T>&);
		matrix<T> & set_row_range(unsigned int,unsigned int,unsigned int,const matrix<T>&);
		matrix<T> & set_chunk(unsigned int,unsigned int,unsigned int,unsigned int,const matrix<T>&);
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
	};
	
	/*Constructors*/
	
	template<class T>
	matrix<T>::matrix(void)
	{
		row = 0;
		col = 0;
		data = nullptr;
	}
	
	template<class T>
	matrix<T>::matrix(unsigned int size)
	{
		if(size < 0) {std::cout << "Runtime Error\nSize cannot be negetive\n"; exit(1);} else {}
		row = size;
		col = size;
		data = new T[size*size];
		if(!data) {std::cout << "Runtime Error\nCannot allocate memory\n"; delete[] data; exit(1);} else {}
	}
	
	template<class T>
	matrix<T>::matrix(unsigned int r,unsigned int c)
	{
		if(r < 0 || c < 0) {std::cout << "Runtime Error\nSize cannot be negetive\n"; exit(1);} else {}
		row = r;
		col = c;
		data = new T[r*c];
		if(!data) {std::cout << "Runtime Error\nCannot allocate memory\n"; delete[] data; exit(1);} else {}
	}
	
	template<class T>
	matrix<T>::matrix(unsigned int size,const T dat[])
	{
		if(size < 0) {std::cout << "Runtime Error\nSize cannot be negetive\n"; exit(1);} else {}
		row = size;
		col = size;
		data = new T[size*size];
		if(!data) {std::cout << "Runtime Error\nCannot allocate memory\n"; delete[] data; exit(1);} else {}

		for(int i = 0; i < size*size; i++)
		{
			data[i] = dat[i];
		}
	}
	
	template<class T>
	matrix<T>::matrix(unsigned int r,unsigned int c,const T dat[])
	{
		if(r < 0 || c < 0) {std::cout << "Runtime Error\nSize cannot be negetive\n"; exit(1);} else {}
		row = r;
		col = c;
		data = new T[r*c];
		if(!data) {std::cout << "Runtime Error\nCannot allocate memory\n"; delete[] data; exit(1);} else {}

		for(int i = 0; i < r*c; i++)
		{
			data[i] = dat[i];
		}
	}

	template<class T>
	matrix<T>::matrix(const matrix &mat)
	{
		row = mat.row;
		col = mat.col;
		data = new T[row*col];
		if(!data) {std::cout << "Runtime Error\nCannot allocate memory\n"; delete[] data; exit(1);} else {}

		for(int i = 0; i < row*col; i++)
		{
			data[i] = mat.data[i];
		}
	}
	
	template<class T>
	matrix<T>::matrix(const T &val)
	{
		row = 1;
		col = 1;
		data = new T[1];
		if(!data) {std::cout << "Runtime Error\nCannot allocate memory\n"; delete[] data; exit(1);} else {}
		data[0] = val;
	}
	
	template<class T>
	matrix<T>::~matrix(void)
	{
		delete[] data;
	}

	/*Printing out the matrix*/
	template<class T>
	std::ostream & operator<<(std::ostream &s,const matrix<T> &mat)
	{
		s << "\n";
		for(int i = 0; i < mat.row; i++)
		{
			s << " ";
			for(int j = 0; j < mat.col; j++)
			{
                s << std::setprecision(5) << std::setw(9) << mat.data[i*mat.col + j] << " ";
			}
			s << "\n";
		}
		s << "\n";
		return s;
	}
	
	/*Indexing*/
	template<class T>
	inline
	unsigned int matrix<T>::if_in(unsigned int size)
	{
		if(size < 0 || size >= row*col) {std::cout << "Runtime Error\nIndex out of bounds\n"; exit(1);} else {}
		return size;
	}
	
	template<class T>
	inline
	unsigned int matrix<T>::if_in(unsigned int r,unsigned int c)
	{
		if (r >= row || c >= col || r < 0 || c < 0) {std::cout << "Runtime Error\nIndex out of bounds\n"; exit(1);} else {}
		return r*col + c;
	}

	template<class T>
	inline
	T & matrix<T>::operator()(unsigned int loc)
	{
		if(row == 1)
		{
			return data[if_in(loc)];
		}
		else if(col == 1)
		{
			return data[if_in(loc)];
		}
		else 
		{
			std::cout << "Runtime Error\nPlease use matrix indexing m(i,j)\n"; exit(1);
		}
	}
	
	template<class T>
	inline
	T & matrix<T>::operator()(unsigned int r,unsigned int c)
	{
		return data[if_in(r,c)];
	}
	
	/*Read only operation for const matrix*/
	
	template<class T>
	inline
	T matrix<T>::operator()(unsigned int loc) const
	{
		if(row == 1)
		{
			if(loc >= row*col) {std::cout << "Runtime Error]\nIndex out of bounds\n"; exit(1);}
			return data[loc];
		}
		else if(col == 1)
		{
			if(loc >= row*col) {std::cout << "Runtime Error]\nIndex out of bounds\n"; exit(1);}
			return data[loc];
		}
		else 
		{
			std::cout << "Runtime Error\nPlease use matrix indexing m(i,j)\n"; exit(1);
		}
	}
	
	template<class T>
	inline
	T matrix<T>::operator()(unsigned int r,unsigned int c) const
	{
		if((r*col + c) >= row*col) {std::cout << "Runtime Error\nIndex out of bounds\n"; exit(1);} else {}
		return data[r*col + c];
	}
	
	/*Read only indexing*/
	template<class T>
	inline
	T matrix<T>::at(unsigned int loc) const
	{
		if(row == 1)
		{
			if(loc >= row*col) {std::cout << "Runtime Error]\nIndex out of bounds\n"; exit(1);}
			return data[loc];
		}
		else if(col == 1)
		{
			if(loc >= row*col) {std::cout << "Runtime Error]\nIndex out of bounds\n"; exit(1);}
			return data[loc];
		}
		else 
		{
			std::cout << "Runtime Error\nPlease use matrix indexing m(i,j)\n"; exit(1);
		}
	}
    
	template<class T>
	inline
	T matrix<T>::at(const unsigned int r,const unsigned int c) const 
	{
		if((r*col + c) >= row*col) {std::cout << "Runtime Error\nIndex out of bounds\n"; exit(1);} else {}
		return data[r*col + c];
	}

	/*Arethemic operations*/
	template<class T>
	matrix<T> & matrix<T>::operator=(const matrix<T> &mat)
	{
		if(row*col != mat.row*mat.col)
		{
			row = mat.row;
			col = mat.col;
			delete[] data;
			data = new T[row*col];
			if(!data) {std::cout << "Runtime Error\nCannot allocate memeory\n"; delete[] data; exit(1);} else {}
			
			for(int i = 0; i < row*col; i++)
			{
				data[i] = mat.data[i];
			}
		}
		else if(row*col == mat.row*mat.col)
		{
			row = mat.row;
			col = mat.col;

			for(int i = 0; i < row*col; i++)
			{
				data[i] = mat.data[i];
			}
		}
		else
		{
			std::cout << "Runtime Error\nUnknown Error occured\n"; delete[] data; exit(1);
		}
		return *this;
	}
	
	template<class T>
	matrix<T> & matrix<T>::operator=(const T &val)
	{
		row = 1;
		col = 1;
		delete[] data;
		data = new T[1];
		if(!data) {std::cout << "Runtime Error\nCannot allocate memeory\n"; delete[] data; exit(1);} else {}
		data[0] = val;
		return *this;
	}
	
	template<class T>
	matrix<T> & matrix<T>::operator+(void)
	{
		return *this;
	}
	
	template<class T>
	matrix<T> & matrix<T>::operator+=(const matrix<T> &mat)
	{
		if((row != mat.row) || (col != mat.col)) {std::cout << "Runtime Error\nMatricies not the same size\n"; exit(1);} else {}

		for(int i = 0; i < row*col; i++)
		{
			data[i] += mat.data[i];
		}
		return *this;
	}
	
	template<class T>
	matrix<T> & matrix<T>::operator-=(const matrix<T> &mat)
	{
		if((row != mat.row) || (col != mat.col)) {std::cout << "Runtime Error\nMatricies not the same size\n"; exit(1);} else {}

		for(int i = 0; i < row*col; i++)
		{
			data[i] -= mat.data[i];
		}
		return *this;
	}

	template<class T>
	matrix<T> & matrix<T>::operator*=(T alpha)
	{
		for(int i = 0; i < row*col; i++)
		{
			data[i] *= alpha;
		}
		return *this;
	}
	
	template<class T>
	matrix<T> & matrix<T>::operator/=(T alpha)
	{
		if(alpha == (T) 0.0) {std::cout << "Runtime Error\nCannot divide by zero\n"; exit(1);} else {}

		for(int i = 0; i < row*col; i++)
		{
			data[i] /= alpha;
		}
		return *this;
	}

	template<class T>
	matrix<T> & matrix<T>::operator&&(const matrix<T> &mat)
	{
		if((row != mat.row) || (col != mat.col)) {std::cout << "Runtime Error\nMatricies not the same size\n"; exit(1);} else {}

		for(int i = 0; i < row*col; i++)
		{
			data[i] *= mat.data[i];
		}
		return *this;
	}
	
	template<class T>
	matrix<T> & matrix<T>::operator||(const matrix<T> &mat)
	{
		if((row != mat.row) || (col != mat.col)) {std::cout << "Runtime Error\nMatricies not the same size\n"; exit(1);} else {}

		for(int i = 0; i < row*col; i++)
		{
			data[i] /= mat.data[i];
		}
		return *this;
	}

	/*Friend arethmeic operations*/
	template<class U>
	matrix<U> operator-(const matrix<U> &mat)
	{
		matrix<U> sol = mat;

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			sol.data[i] *= (U) -1.0;
		}
		return sol;
	}
	
	template<class U>
	matrix<U> operator-(const matrix<U> &mat1,const matrix<U> &mat2)
	{
		if(mat1.row != mat2.row || mat1.col != mat2.col) {std::cout << "Runtime Error\nMatricies not the same size\n"; exit(1);} else {}
		matrix<U> sol(mat1.row,mat1.col);

		for(int i = 0; i < mat1.row*mat2.col; i++)
		{
			sol.data[i] = mat1.data[i] - mat2.data[i];
		}
		return sol;
	}
	
	template<class U>
	matrix<U> operator+(const matrix<U> &mat1,const matrix<U> &mat2)
	{
		if(mat1.row != mat2.row || mat1.col != mat2.col) {std::cout << "Runtime Error\nMatricies not the same size\n"; exit(1);} else {}
		matrix<U> sol(mat1.row,mat1.col);

		for(int i = 0; i < mat1.row*mat1.col; i++)
		{
			sol.data[i] = mat1.data[i] + mat2.data[i];
		}
		return sol;
	}

	template<class U>
	matrix<U> operator*(U alpha,const matrix<U> &mat)
	{
		matrix<U> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			sol.data[i] = alpha*mat.data[i];
		}
		return sol;
	}

	template<class U>
	matrix<U> operator*(const matrix<U> &mat,U alpha)
	{
		matrix<U> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			sol.data[i] = alpha*mat.data[i];
		}
		return sol;
	}
	
	template<class U>
	matrix<U> operator/(U alpha,const matrix<U> &mat)
	{
		matrix<U> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			sol.data[i] = alpha/mat.data[i];
		}
		return sol;
	}

	template<class U>
	matrix<U> operator/(const matrix<U> &mat,U alpha)
	{
		matrix<U> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			sol.data[i] = mat.data[i]/alpha;
		}
		return sol;
	}

	template<class U>
	matrix<U> operator&(const matrix<U> &mat1,const matrix<U> &mat2)
	{
		if(mat1.row != mat2.row || mat1.col != mat2.col) {std::cout << "Runtime Error\nMatricies not the same size\n"; exit(1);} else {}
		matrix<U> sol(mat1.row,mat2.col);

		for(int i = 0; i < mat1.row*mat2.col; i++)
		{
			sol.data[i] = mat1.data[i]*mat2.data[i];
		}
		return sol;
	}

	template<class U>
	matrix<U> operator|(const matrix<U> &mat1,const matrix<U> &mat2)
	{
		if(mat1.row != mat2.row || mat1.col != mat2.col) {std::cout << "Runtime Error\nMatricies not the same size\n"; exit(1);} else {}
		matrix<U> sol(mat1.row,mat2.col);

		for(int i = 0; i < mat1.row*mat2.col; i++)
		{
			sol.data[i] = mat1.data[i]/mat2.data[i];
		}
		return sol;
	}

	/*Friend generating functions*/
	template<class U>
	matrix<U> matrix<U>::zeros(unsigned int size)
	{
		matrix<U> sol(size);

		for(int i = 0; i < size*size; i++)
		{
			sol.data[i] = (U) 0.0;
		}
		return sol;
	}
	
	template<class U>
	matrix<U> matrix<U>::zeros(unsigned int r,unsigned int c)
	{
		matrix<U> sol(r,c);

		for(int i = 0; i < r*c; i++)
		{
			sol.data[i] = (U) 0.0;
		}
		return sol;
	}
	
	template<class U>
	matrix<U> matrix<U>::ones(unsigned int size)
	{
		matrix<U> sol(size);

		for(int i = 0; i < size*size; i++)
		{
			sol.data[i] = (U) 1.0;
		}
		return sol;
	}
	
	template<class U>
	matrix<U> matrix<U>::ones(unsigned int r,unsigned int c)
	{
		matrix<U> sol(r,c);

		for(int i = 0; i < r*c; i++)
		{
			sol.data[i] = (U) 1.0;
		}
		return sol;
	}
	
	/*Uninitalized matrix*/
	template<class T>
	matrix<T> matrix<T>::shape(unsigned int size)
	{
		matrix<T> sol(size);
		return sol;
	}
	
	template<class T>
	matrix<T> matrix<T>::shape(unsigned int r,unsigned int c)
	{
		matrix<T> sol(r,c);
		return sol;
	}
	
	template<class U>
	matrix<U> matrix<U>::randu(unsigned int size)
	{
		matrix<U> sol(size);

		for(int i = 0; i < size*size; i++)
		{
			sol.data[i] = (U) (double)rand()/(double)RAND_MAX;
		}
		return sol;
	}
	
	template<class U>
	matrix<U> matrix<U>::randu(unsigned int r,unsigned int c)
	{
		matrix<U> sol(r,c);

		for(int i = 0; i < r*c; i++)
		{
			sol.data[i] = (U) (double)rand()/(double)RAND_MAX;
		}
		return sol;
	}
	
	template<class U>
	matrix<U> matrix<U>::randu(double a,double b,unsigned int size)
	{
		matrix<U> sol(size);

		for(int i = 0; i < size*size; i++)
		{
			sol.data[i] = (U) ((b-a)*(double)rand()/(double)RAND_MAX + a);
		}
		return sol;
	}
	
	template<class U>
	matrix<U> matrix<U>::randu(double a,double b,unsigned int r,unsigned int c)
	{
		matrix<U> sol(r,c);

		for(int i = 0; i < r*c; i++)
		{
			sol.data[i] = (U) ((b-a)*(double)rand()/(double)RAND_MAX + a);
		}
		return sol;
	}
	
	template<class U>
	matrix<U> matrix<U>::eye(unsigned int size)
	{
		matrix<U> sol = matrix<U>::zeros(size);
 
		for(int i = 0; i < size; i++)
		{
			sol.data[i*size + i] = (U) 1.0;
		}
		return sol;
	}
	
	template<class U>
	matrix<U> matrix<U>::linspace(U a,U b,unsigned int size)
	{
		if(a > b) {U temp = a; a = b; b = temp;} /*Swaps if range is negetive*/
		else if(b > a) {/*Do nothing*/}
		else {std::cout << "Runtime Error\na & b cannot be equal\n"; exit(1);}
		U h = (U) (b-a)/(U)(size-1);
		matrix<U> sol(1,size);
		sol(0) = (U) a;
		sol(size-1) = (U) b;
		for(int i = 1; i < size-1; i++)
		{
			sol(i) = (U) (sol(i-1) + h);
		}
		return sol;
	}
	
	template<class U>
	matrix<U> matrix<U>::delspace(U a,U b,U h)
	{
		if(a > b) {U temp = a; a = b; b = temp;} /*Swaps if range is negetive*/
		else if(b > a) {/*Do nothing*/}
		else {std::cout << "Runtime Error\na & b cannot be equal\n"; exit(1);}
		//double h = (double) (b-a)/(size-1);
		unsigned int size = (b-a)/h + (U)1; 
		matrix<U> sol(1,size);
		sol(0) = (U) a;
		sol(size-1) = (U) b;
		for(int i = 1; i < size-1; i++)
		{
			sol(i) = (U) (sol(i-1) + h);
		}
		return sol;
	}
	
	/*Normal random variables*/
	template<class T>
	matrix<T> matrix<T>::randn(unsigned int N)
	{
		matrix<T> sol(N);

		for(int i = 0; i < sol.row*sol.col; ++i)
		{
			double V = (double) (((double)2)*(double)rand()/(double)RAND_MAX - (double)1);
			double U = (double) (((double)2)*(double)rand()/(double)RAND_MAX - (double)1);
			double S = (double) std::pow(U,2) + std::pow(V,2);
		
			while(S >= (double)1)
			{
				V = (double) (((double)2)*(double)rand()/(double)RAND_MAX - (double)1);
				U = (double) (((double)2)*(double)rand()/(double)RAND_MAX - (double)1);
				S = (double) std::pow(U,2) + std::pow(V,2);
			}
			sol.data[i] = (T) U*std::sqrt(-(double)2*std::log(S)/S);
		}
		return sol;
	}
	
	template<class T>
	matrix<T> matrix<T>::randn(unsigned int r,unsigned int c)
	{
		matrix<T> sol(r,c);

		for(int i = 0; i < sol.row*sol.col; ++i)
		{
			double V = (double) (((double)2)*(double)rand()/(double)RAND_MAX - (double)1);
			double U = (double) (((double)2)*(double)rand()/(double)RAND_MAX - (double)1);
			double S = (double) std::pow(U,2) + std::pow(V,2);
		
			while(S >= (double)1)
			{
				V = (double) (((double)2)*(double)rand()/(double)RAND_MAX - (double)1);
				U = (double) (((double)2)*(double)rand()/(double)RAND_MAX - (double)1);
				S = (double) std::pow(U,2) + std::pow(V,2);
			}
			sol.data[i] = (T) U*std::sqrt(-(double)2*std::log(S)/S);
		}
		return sol;
	}

	template<class T>
	matrix<T> matrix<T>::randn(double mu,double var,unsigned int N)
	{
		matrix<T> sol(N);

		for(int i = 0; i < sol.row*sol.col; ++i)
		{
			double V = (double) (((double)2)*(double)rand()/(double)RAND_MAX - (double)1);
			double U = (double) (((double)2)*(double)rand()/(double)RAND_MAX - (double)1);
			double S = (double) std::pow(U,2) + std::pow(V,2);
		
			while(S >= (double)1)
			{
				V = (double) (((double)2)*(double)rand()/(double)RAND_MAX - (double)1);
				U = (double) (((double)2)*(double)rand()/(double)RAND_MAX - (double)1);
				S = (double) std::pow(U,2) + std::pow(V,2);
			}
			sol.data[i] = (T) (mu + std::sqrt(var)*(U*std::sqrt(-(double)2*std::log(S)/S)));
		}
		return sol;
	}
	
	template<class T>
	matrix<T> matrix<T>::randn(double mu,double var,unsigned int r,unsigned int c)
	{
		matrix<T> sol(r,c);

		for(int i = 0; i < sol.row*sol.col; ++i)
		{
			double V = (double) (((double)2)*(double)rand()/(double)RAND_MAX - (double)1);
			double U = (double) (((double)2)*(double)rand()/(double)RAND_MAX - (double)1);
			double S = (double) std::pow(U,2) + std::pow(V,2);
		
			while(S >= (double)1)
			{
				V = (double) (((double)2)*(double)rand()/(double)RAND_MAX - (double)1);
				U = (double) (((double)2)*(double)rand()/(double)RAND_MAX - (double)1);
				S = (double) std::pow(U,2) + std::pow(V,2);
			}
			sol.data[i] = (T) (mu + std::sqrt(var)*(U*std::sqrt(-(double)2*std::log(S)/S)));
		}
		return sol;
	}
	template<class T>
	void matrix<T>::set_rand(void)
	{
		srand(time(NULL)*123456789);
	}
	
	template<class T>
	void matrix<T>::set_rand(long long val)
	{
		srand(val);
	}
	
	/*Matrix - vector operations*/
	template<class T>
	matrix<T> matrix<T>::t(void)
	{
		matrix<T> trans(col,row);

		for(int i = 0; i < row; i++)
		{
			for(int j = 0; j < col; j++)
			{
				trans.data[j*row + i] = data[i*col + j];
			}
		}
		return trans;
	}
	
	/*Const transpose - used for const mats*/
	template<class T>
	matrix<T> matrix<T>::t(void) const
	{
		matrix<T> trans(col,row);

		for(int i = 0; i < row; i++)
		{
			for(int j = 0; j < col; j++)
			{
				trans.data[j*row + i] = data[i*col + j];
			}
		}
		return trans;
	}
	
	template<class U>
	matrix<U> row_concat(const matrix<U> &mat1,const matrix<U> &mat2)
	{
		if(mat1.col != mat2.col) {std::cout << "Runtime Error\nUnequal columns\n"; exit(1);} else {}
		matrix<U> sol(mat1.row+mat2.row,mat1.col);
		
		for(int i = 0; i < mat1.row*mat1.col; i++)
		{
			sol.data[i] = mat1.data[i];
		}
		
		int position = mat1.row*mat1.col;

		for(int i = 0; i < mat2.row*mat2.col; i++)
		{
			sol.data[i + position] = mat2.data[i];
		}
		return sol;
	}
	
	template<class U>
	matrix<U> col_concat(const matrix<U> &m1,const matrix<U> &m2)
	{
		matrix<U> m1_cpy = m1;
		matrix<U> m2_cpy = m2;
		return row_concat(m1_cpy.t(),m2_cpy.t()).t();
	}
	
	template<class T>
	inline
	unsigned int matrix<T>::rows(void) const
	{
		return row;
	}
	
	template<class T>
	inline
	unsigned int matrix<T>::cols(void) const
	{
		return col;
	}
	
	template<class T>
	inline
	unsigned int matrix<T>::numel(void) const 
	{
		return col*row;
	}
	
	template<class T>
	matrix<T> matrix<T>::vectorize(void)
	{
		matrix<T> sol(row*col,1);
		sol.row = row*col;
		sol.col = 1;

		for(int i = 0; i < row*col; i++)
		{
			sol.data[i] = data[i];
		}
		return sol;
	}
	
	template<class T>
	matrix<T> & matrix<T>::fill(T alpha)
	{
		for(int i = 0; i < row*col; i++)
		{
			data[i] = alpha;
		}
		return *this;
	}

	template<class T>
	matrix<T> matrix<T>::vectorize(unsigned int c)
	{
		if(c == 0)
		{
			matrix<T> sol(row*col,1);
			sol.row = row*col;
			sol.col = 1;

			for(int i = 0; i < row*col; i++)
			{
				sol.data[i] = data[i];
			}
			return sol;
		}
		else if(c == 1)
		{
			matrix<T> sol(1,row*col);
			sol.col = row*col;
			sol.row = 1;

			for(int i = 0; i < row*col; i++)
			{
				sol.data[i] = data[i];
			}
			return sol;
		}
		else 
		{
			std::cout << "Runtime Error\nVectorize only supports values 0 or 1\n"; exit(1);
		}
	}
	
	template<class T>
	matrix<T> & matrix<T>::for_row(unsigned int r,T (*f)(T))
	{
		if(r >= row || r < 0) {std::cout << "Runtime Error\nIndex out of bounds\n"; exit(1);} else {}

		for(int j = 0; j < col; j++)
		{
			data[r*col + j] = f(data[r*col + j]);
		}
		return *this;	
	}
	
	template<class T>
	matrix<T> & matrix<T>::for_col(unsigned int c,T (*f)(T))
	{
		if(c >= col || c < 0) {std::cout << "Runtime Error\nIndex out of bounds\n"; exit(1);} else {}

		for(int i = 0; i < row; i++)
		{
			data[i*col + c] = f(data[i*col + c]);
		}
		return *this;
	}
	
	template<class T>
	matrix<T> & matrix<T>::for_col_range(unsigned int c,unsigned int a,unsigned int b,T (*f)(T))
	{
		if(c >= col || c < 0 || a >= row || a < 0 || b >= row || b < 0) {std::cout << "Runtime Error\nIndex out of range\n"; exit(1);}
		if(b < a) {unsigned int temp = a; a = b; a = b; b = temp;} else {/*Do nothing*/}
		
		for(int i = a; i <=b; i++)
		{
			data[i*col + c] = f(data[i*col + c]);
		}
		return *this;
	}
	
	template<class T>
	matrix<T> & matrix<T>::for_row_range(unsigned int r,unsigned int a,unsigned int b,T (*f)(T))
	{
		if(r >= row || r < 0 || a >= col || a < 0 || b >= col || b < 0) {std::cout << "Runtime Error\nIndex out of range\n"; exit(1);}
		if(b < a) {unsigned int temp = a; a = b; a = b; b = temp;} else {/*Do nothing*/}
		
#pragma omp parallel for
		for(int i = a; i<=b; i++)
		{
			data[r*col + i] = f(data[r*col + i]);
		}
		return *this;
	}
	
	template<class T>
	matrix<T> & matrix<T>::for_all(T (*f)(T))
	{
		for(int i = 0; i < row*col; i++)
		{
			data[i] = f(data[i]);
		}
		return *this;
	}
	
	template<class T>
	matrix<T> & matrix<T>::row_swap(unsigned int r1,unsigned int r2)
	{
		if(r1 >= row || r1 < 0 || r2 >= row || r2 < 0) {std::cout << "Runtime Error\nIndex out of bounds\n"; exit(1);} else {}
		T *temp = new T[col];
		T *temp2 = new T[col];
		if(!temp || !temp2) {std::cout << "Runtime Error\nCannot allocate memory\n";delete[] temp,temp2; exit(1);} else {}
		
		for(int j = 0; j < col; j++)
		{
			temp[j] = data[r1*col + j];
			temp2[j] = data[r2*col + j];
		}

		for(int j = 0; j < col; j++)
		{
			data[r1*col + j] = temp2[j];
			data[r2*col + j] = temp[j];
		}
		delete[] temp2,temp;
		return *this;
	}
	
	template<class T>
	matrix<T> & matrix<T>::col_swap(unsigned int c1,unsigned int c2)
	{
		if(c1 >= col || c1 < 0 || c2 >= col || c2 < 0) {std::cout << "Runtime Error\nIndex out of bounds\n"; exit(1);} else {}
		T *temp = new T[row];
		T *temp2 = new T[row];
		if(!temp || !temp2) {std::cout << "Runtime Error\nCannot allocate memory\n";delete[] temp,temp2; exit(1);} else {}
		
		for(int i = 0; i < row; i++)
		{
			temp[i] = data[i*col + c1];
			temp2[i] = data[i*col + c2];
		}

		for(int i = 0; i < row; i++)
		{
			data[i*col + c1] = temp2[i];
			data[i*col + c2] = temp[i];
		}
		delete[] temp2,temp;
		return *this;
	}
	
	/*Returns the row or col as a new vector*/
	template<class T>
	matrix<T> matrix<T>::get_row(unsigned int r)
	{
		if(r >= row || r < 0) {std::cout << "Runtime Error\nIndex out of bounds\n"; exit(1);} else {}
		matrix<T> sol(1,col);
		
		for(int i = 0; i < col; i++)
		{
			sol.data[i] = data[r*col + i];
		}
		return sol;
	}
	
	template<class T>
	matrix<T> matrix<T>::get_col(unsigned int c)
	{
		if(c >= col || c < 0) {std::cout << "Runtime Error\nIndex out of bounds\n"; exit(1);} else {}
		matrix<T> sol(row,1);
		
		for(int i = 0; i < row; i++)
		{
			sol.data[i] = data[i*col + c];
		}
		return sol;
	}
	
	template<class T>
	matrix<T> matrix<T>::get_row_range(unsigned int r,unsigned int a,unsigned int b)
	{
		if(r >= row || r < 0 || a >= col || a < 0 || b >= col || b < 0) {std::cout << "Runtime Error\nIndex out of bounds\n"; exit(1);} else {}
		if(b < a) {unsigned int temp = b; b = a; a = temp;} else {} /*Swaps the range if it is invalid*/
		matrix<T> sol(1,b-a + 1);
		int j = 0;
#pragma omp for
		for(int i = a; i <= b; i++)
		{
			sol.data[j++] = data[r*col + i];
		}
		return sol;
	}
	
	template<class T>
	matrix<T> matrix<T>::get_col_range(unsigned int c,unsigned int a,unsigned int b)
	{
		if(c >= col || c < 0 || a >= row || a < 0 || b >= row || b < 0) {std::cout << "Runtime Error\nIndex out of bounds\n";exit(1);} else {}
		if(b < a) {unsigned int temp = b; b = a; a = temp;} else {} /*Swaps the range if it is invalid*/
		matrix<T> sol(b-a+1,1);
		int j = 0;
#pragma omp for
		for(int i = a; i <= b; i++)
		{
			sol.data[j++] = data[i*col + c];
		}
		return sol;
	}
	
	template<class T>
	matrix<T> matrix<T>::diag(void)
	{
		if(row != col) {std::cout << "Runtime Error\nMatrix must be square\n"; exit(1);} else {}
		matrix<T> d(row,1);

		for(int i = 0; i < row; i++)
		{
			d.data[i] = data[i*col + i];
		}
		return d;
	}
	
	template<class T>
	unsigned int matrix<T>::is_in(T alpha)
	{
		for(int i = 0; i < row*col; i++)
		{
			if(data[i] == alpha) {return 1;}
			else {continue;}
		}
		return 0;
	}
	
	/*Boolean tests for the shape*/
	template<class T>
	inline
	unsigned int matrix<T>::is_square(void) const
	{
		if(row == col) {return 1;}
		else if (row != col) {return 0;}
		else {std::cout << "Runtime Error\nSomething went wrong\n"; exit(1);}
	}
	
	template<class T>
	inline
	unsigned int matrix<T>::is_vector(void) const
	{
		if(row == 1 || col == 1) {return 1;}
		else if(row != 1 || col != 1) {return 0;}
		else {std::cout << "Runtime Error\nSomething went wrong\n"; exit(1);}
	}
	
	template<class T>
	inline
	unsigned int matrix<T>::is_row_vector(void) const
	{
		if(row == 1) {return 1;}
		else if(row != 1) {return 0;}
		else {std::cout << "Runtime Error\nSomething went wrong\n"; exit(1);}
	}
	
	template<class T>
	inline
	unsigned int matrix<T>::is_col_vector(void) const
	{
		if(col == 1) {return 1;}
		else if(col != 1) {return 0;}
		else {std::cout << "Runtime Error\nSomething went wrong\n"; exit(1);}
	}
	
	/*Matrix Multiplication*/
	template<class U>
	matrix<U> operator*(const matrix<U> &mat1,const matrix<U> &mat2)
	{	
		if(mat1.col != mat2.row) {std::cout << "Runtime Error\nMatrix not equivlent\n"; exit(1);} else {}
		matrix<U> sol = matrix<U>::zeros(mat1.rows(),mat2.cols());
#pragma omp parallel for
		for(int i = 0; i < mat1.row; i++)
		{
			for(int k = 0; k < mat1.col; k++)
			{
				U a = mat1.data[i*mat1.col + k];
				for(int j = 0; j < mat2.col; j++)
				{
					sol.data[i*sol.col + j] += a*mat2.data[k*mat2.col + j];
				}
			}
		}
		return sol;
	}
	
	/*Uses multiplication for banded matricies*/
	template<class U>
	matrix<U> mat_mull_vec_band(const matrix<U> &mat1,const matrix<U> &mat2,unsigned int band_len)
	{
		if(mat1.col != mat2.row || !mat2.is_col_vector()) {std::cout << "Runtime Error\nMatrix not the correct size or mat2 not col vector\n"; exit(1);} else {}
		matrix<U> sol = matrix<U>::zeros(mat1.row,1);
		
		/*Solve the sparse system*/
		/*Internal part of matrix*/
#pragma omp parallel for
		for(int i = band_len; i < mat1.row-band_len; i++)
		{
			for(int k = i-band_len; k < i+band_len; k++)
			{
				sol.data[i] += mat1.data[i*mat1.col + k]*mat2.data[k];
			}
		}
		
		/*Top of the matrix*/
#pragma omp parallel for
		for(int i = 0; i < band_len; i++)
		{
			for(int k = 0; k < i+band_len; k++)
			{
				sol.data[i] += mat1.data[i*mat1.col + k]*mat2.data[k];
			}
		}
		
		/*bottom part of the matrix*/
#pragma omp parallel for
		for(int i = mat1.row-band_len; i < mat1.row; i++)
		{
			for(int k = i-band_len; k < mat1.row; k++)
			{
				sol.data[i] += mat1.data[i*mat1.col + k]*mat2.data[k];
			}
		}
		return sol;
	}
	
	/*Stats functions! -> more to come later when needed!*/
	template<class U>
	U sum(const matrix<U> &mat)
	{
		U s = (U) 0.0;

		for(int i = 0; i < mat.row*mat.col; i++)
			{
				s += mat.data[i];
			}
			return s;
	}
	
	template<class U>
	U mean(const matrix<U> &mat)
	{
		return sum(mat)/(U)(mat.row*mat.col);
	}
	
	template<class U>
	U var(const matrix<U> &mat)
	{
		U v = (U) 0.0;
		U m = mean(mat);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			v += std::pow(mat.data[i] - m,(U)2);
		}
		return v/((U)mat.row*mat.col - 1);
	}
	
	template<class U>
	U sd(const matrix<U> &mat)
	{
		return std::sqrt(var(mat));
	}
	
	template<class U>
	matrix<U> cov(const matrix<U> &A,const matrix<U> &B)
	{
		if((A.numel() != B.numel()) || (!A.is_vector()) || (!B.is_vector())) {std::cout << "Runtime Error\nMatricies not the correct size\n"; exit(1);} else {}
		U A_mean = mean(A);
		U B_mean = mean(B);
		U A_var = var(A);
		U B_var = var(B);
		matrix<U> sol = matrix<U>::zeros(2,2);
		
		U c = (U) 0;

		for(int i = 0; i < A.numel(); i++)
		{
			c += (A.data[i] - A_mean)*(B.data[i] - B_mean);
		}
		c /= (U)(A.numel() - 1);
		
		sol(0,0) = A_var;
		sol(0,1) = c;
		sol(1,0) = c;
		sol(1,1) = B_var;
		return sol;
	}
	
	template<class U>
	matrix<U> cov(const matrix<U> &X)
	{
		/*Initalize our covariance matrix*/
		matrix<U> C = matrix<U>::zeros(X.rows(),X.rows());
		/*Calculate the means ahead of time*/
		matrix<U> M = matrix<U>::zeros(X.rows(),1); 
#pragma omp parallel for
		for(int i = 0; i < X.rows(); i++)
		{
			U m = (U)0;
			for(int j = 0; j < X.cols(); j++)
			{
				m += X(i,j);
			}
			M(i) = m/(U)X.cols();
		}
#pragma omp parallel for 
		for(int i = 0; i < X.rows(); i++)
		{
			for(int j = 0; j < X.rows(); j++)
			{
				for(int k = 0; k < X.cols(); k++)
				{
					C(i,j) += (X(i,k) - M(i))*(X(j,k) - M(j));
				}
				C(i,j) = C(i,j)/((U)X.cols() - (U)1);
			}
		}
		return C;
	}
	
	template<class U>
	U max(const matrix<U> &mat)
	{
		U m = mat.data[0];
		
		for(int i = 1; i < mat.row*mat.col; i++)
		{
			if(mat.data[i] > m) {m = mat.data[i];} else {}
		}
		return m;
	}
	
	template<class U>
	U min(const matrix<U> &mat)
	{
		U m = mat.data[0] /*Large constant*/;

		for(int i = 1; i < mat.row*mat.col; i++)
		{
			if(mat.data[i] < m) {m = mat.data[i];} else {}
		}
		return m;
	}
	
	template<class T>
	matrix<T> matrix<T>::randperm(unsigned int N)
	{
		matrix<T> gen = matrix<T>::zeros(1,N);
#pragma omp for
		for(int i = 0; i < N; i++)
		{
			T rand_ = (T) (((double)N+1.)*(double)rand()/(double)RAND_MAX + ((double)N+1.)*(double)rand()/(double)RAND_MAX)/(double)2;
			while(sum(gen == rand_) >= 1.)
			{
				rand_ = (T) (((double)N+1.)*(double)rand()/(double)RAND_MAX + ((double)N+1.)*(double)rand()/(double)RAND_MAX)/(double)2;
			}
			gen.data[i] = rand_;
		}
		/*Does random col swaps to add more randomness*/
#pragma omp for
		for(int i = 0; i < N; i++)
		{
			T swap_coe_1 = (T) (((double)N)*(double)rand()/(double)RAND_MAX + ((double)N)*(double)rand()/(double)RAND_MAX)/(double)2;
			T swap_coe_2 = (T) (((double)N)*(double)rand()/(double)RAND_MAX + ((double)N)*(double)rand()/(double)RAND_MAX)/(double)2;
			gen.col_swap(swap_coe_1,swap_coe_2);
		}
		return gen;
	}
	
	template<class T>
	matrix<T> matrix<T>::randperm(unsigned int N,unsigned int size)
	{
		if(size > N) {std::cout << "Runtime Error\nSize cannot be greater then the length of the vector\n"; exit(1);} else {}
		matrix<T> gen = matrix<T>::randperm(N);
		matrix<T> sol(1,size);
	
		/*Does random col swaps to add more randomness*/
#pragma omp for
		for(int i = 0; i < N; i++)
		{
			T swap_coe_1 = (T) (((double)N)*(double)rand()/(double)RAND_MAX + ((double)N)*(double)rand()/(double)RAND_MAX)/(double)2;
			T swap_coe_2 = (T) (((double)N)*(double)rand()/(double)RAND_MAX + ((double)N)*(double)rand()/(double)RAND_MAX)/(double)2;
			gen.col_swap(swap_coe_1,swap_coe_2);
		}
		/*Get the size of what we want to return*/

		for(int i = 0; i < size; i++)
		{
			sol.data[i] = gen.data[i];
		}
		return sol;
	}
	
	/*More stat functions to come later*/
	
	/*General functions*/
	template<class U>
	U inf_norm(const matrix<U> &mat)
	{
		matrix<U> mat_cpy = mat;
		U s = (U) 0.0;
		U t = (U) 0.0;
		for(int i = 0; i < mat_cpy.row; i++)
		{
			s = sum(mat_cpy.get_row(i));
			if(s > t) {t = s;}
		}
		return t;
	}
	
	template<class U>
	U one_norm(const matrix<U> &mat)
	{
		matrix<U> mat_cpy = mat;
		return inf_norm(mat_cpy.t());
	}
	
	template<class U>
	U norm_2(const matrix<U> &mat)
	{
		/*If a matrix just creates a bound*/
		U n = (U) 0.0;

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			n += std::pow(mat.data[i],(double)2.0);
		}
		return (U) std::sqrt(n);
	}
	
	template<class U>
	U p_norm(const matrix<U> &mat,double alpha)
	{
		U n = (U) 0.0;

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			n += std::pow(mat.data[i],alpha);
		}
		return std::pow(n,(double)1.0/alpha);
	}
	
	template<class U>
	U rcond(const matrix<U> &mat)
	{
		matrix<U> mat_cpy = mat;
		matrix<U> temp = mat_cpy.diag();
		return (U) 1.0/(max(temp)/min(temp));
	}
	
	/*Cmath functions*/
	template<class U>
	matrix<U> cos(const matrix<U> &mat)
	{
		matrix<U> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			sol.data[i] = std::cos(mat.data[i]);
		}
		return sol;
	}
	
	template<class U>
	matrix<U> sin(const matrix<U> &mat)
	{
		matrix<U> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			sol.data[i] = std::sin(mat.data[i]);
		}
		return sol;
	}
	
	template<class U>
	matrix<U> tan(const matrix<U> &mat)
	{
		matrix<U> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			sol.data[i] = std::tan(mat.data[i]);
		}
		return sol;
	}
	
	template<class U>
	matrix<U> acos(const matrix<U> &mat)
	{
		matrix<U> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			sol.data[i] = std::acos(mat.data[i]);
		}
		return sol;
	}
	
	template<class U>
	matrix<U> asin(const matrix<U> &mat)
	{
		matrix<U> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			sol.data[i] = std::asin(mat.data[i]);
		}
		return sol;
	}
	
	template<class U>
	matrix<U> atan(const matrix<U> &mat)
	{
		matrix<U> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			sol.data[i] = std::atan(mat.data[i]);
		}
		return sol;
	}
	
	template<class U>
	matrix<U> cosh(const matrix<U> &mat)
	{
		matrix<U> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			sol.data[i] = std::cosh(mat.data[i]);
		}
		return sol;
	}
	
	template<class U>
	matrix<U> sinh(const matrix<U> &mat)
	{
		matrix<U> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			sol.data[i] = std::sinh(mat.data[i]);
		}
		return sol;
	}
	
	template<class U>
	matrix<U> tanh(const matrix<U> &mat)
	{
		matrix<U> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			sol.data[i] = std::tanh(mat.data[i]);
		}
		return sol;
	}
	
	template<class U>
	matrix<U> acosh(const matrix<U> &mat)
	{
		matrix<U> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			sol.data[i] = std::acosh(mat.data[i]);
		}
		return sol;
	}
	
	template<class U>
	matrix<U> asinh(const matrix<U> &mat)
	{
		matrix<U> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			sol.data[i] = std::asinh(mat.data[i]);
		}
		return sol;
	}
	
	template<class U>
	matrix<U> atanh(const matrix<U> &mat)
	{
		matrix<U> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			sol.data[i] = std::atanh(mat.data[i]);
		}
		return sol;
	}
	
	template<class U>
	matrix<U> exp(const matrix<U> &mat)
	{
		matrix<U> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			sol.data[i] = std::exp(mat.data[i]);
		}
		return sol;
	}
	
	template<class U>
	matrix<U> log(const matrix<U> &mat)
	{
		matrix<U> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			sol.data[i] = std::log(mat.data[i]);
		}
		return sol;
	}
	
	template<class U>
	matrix<U> log10(const matrix<U> &mat)
	{
		matrix<U> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			sol.data[i] = std::log10(mat.data[i]);
		}
		return sol;
	}
	
	template<class U>
	matrix<U> log2(const matrix<U> &mat)
	{
		matrix<U> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			sol.data[i] = std::log2(mat.data[i]);
		}
		return sol;
	}
	
	/*Vectorized power method for convience*/
	template<class U>
	matrix<U> operator^(const matrix<U> &mat,double alpha)
	{
		matrix<U> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			sol.data[i] = std::pow(mat.data[i],alpha);
		}
		return sol;
	}
	
	template<class U>
	matrix<U> operator^(double alpha,const matrix<U> &mat)
	{
		matrix<U> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			sol.data[i] = std::pow(alpha,mat.data[i]);
		}
		return sol;
	}
	
	template<class U>
	matrix<U> ceil(const matrix<U> &mat)
	{
		matrix<U> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			sol.data[i] = std::ceil(mat.data[i]);
		}
		return sol;
	}
	
	template<class U>
	matrix<U> floor(const matrix<U> &mat)
	{
		matrix<U> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			sol.data[i] = std::floor(mat.data[i]);
		}
		return sol;
	}
	
	template<class U>
	matrix<U> abs(const matrix<U> &mat)
	{
		matrix<U> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			sol.data[i] = std::abs(mat.data[i]);
		}
		return sol;
	}
	
	/*Calulates the diff in column wise ordering (just transpose if you want different)*/
	template<class U>
	matrix<U> diff(const matrix<U> &mat)
	{
		matrix<U> d(mat.row,mat.col-1);

		for(int i = 1; i < mat.row*(mat.col); i++)
		{
			d.data[i-1] = mat.data[i] - mat.data[i-1];
		}
		return d;
	}
	
	/*Similar matlab operators for convience*/
	
	template<class U>
	matrix<U> operator+(const matrix<U> &mat,U alpha)
	{
		matrix<U> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			sol.data[i] = mat.data[i] + alpha;
		}
		return sol;
	}
	
	template<class U>
	matrix<U> operator+(U alpha,const matrix<U> &mat)
	{
		matrix<U> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			sol.data[i] = mat.data[i] + alpha;
		}
		return sol;
	}
	
	template<class U>
	matrix<U> operator-(U alpha,const matrix<U> &mat)
	{
		matrix<U> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			sol.data[i] = alpha - mat.data[i];
		}
		return sol;
	}
	
	template<class U>
	matrix<U> operator-(const matrix<U> &mat,U alpha)
	{
		matrix<U> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			sol.data[i] = mat.data[i] - alpha;
		}
		return sol;
	}
	
	template<class T>
	matrix<T> & matrix<T>::operator++(void)
	{
		for(int i = 0; i < row*col; i++)
		{
			data[i] += (T)1;
		}
		return *this;
	}
	/*Logical operations on a mtrix*/
	
	template<class U>
	matrix<unsigned int> operator>(const matrix<U> &mat,U alpha)
	{
		matrix<unsigned int> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			if(mat.data[i] > alpha) {sol.data[i] = 1;}
			else if(mat.data[i] <= alpha) {sol.data[i] = 0;}
			else{std::cout << "Runtime Error\nSomething went wrong with logic\n"; exit(1);}
		}
		return sol;
	}
	
	template<class U>
	matrix<unsigned int> operator>(U alpha,const matrix<U> &mat)
	{
		matrix<unsigned int> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			if(alpha > mat.data[i]) {sol.data[i] = 1;}
			else if(alpha <= mat.data[i]) {sol.data[i] = 0;}
			else{std::cout << "Runtime Error\nSomething went wrong with logic\n"; exit(1);}
		}
		return sol;
	}
	
	template<class U>
	matrix<unsigned int> operator<(const matrix<U> &mat,U alpha)
	{
		matrix<unsigned int> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			if(mat.data[i] < alpha) {sol.data[i] = 1;}
			else if(mat.data[i] >= alpha) {sol.data[i] = 0;}
			else{std::cout << "Runtime Error\nSomething went wrong with logic\n"; exit(1);}
		}
		return sol;
	}
	
	template<class U>
	matrix<unsigned int> operator<(U alpha,const matrix<U> &mat)
	{
		matrix<unsigned int> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			if(alpha < mat.data[i]) {sol.data[i] = 1;}
			else if(alpha >= mat.data[i]) {sol.data[i] = 0;}
			else{std::cout << "Runtime Error\nSomething went wrong with logic\n"; exit(1);}
		}
		return sol;
	}
	
	/*Equal operators*/
	template<class U>
	matrix<unsigned int> operator>=(const matrix<U> &mat,U alpha)
	{
		matrix<unsigned int> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			if(mat.data[i] >= alpha) {sol.data[i] = 1;}
			else if(mat.data[i] < alpha) {sol.data[i] = 0;}
			else{std::cout << "Runtime Error\nSomething went wrong with logic\n"; exit(1);}
		}
		return sol;
	}
	
	template<class U>
	matrix<unsigned int> operator>=(U alpha,const matrix<U> &mat)
	{
		matrix<unsigned int> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			if(alpha >= mat.data[i]) {sol.data[i] = 1;}
			else if(alpha < mat.data[i]) {sol.data[i] = 0;}
			else{std::cout << "Runtime Error\nSomething went wrong with logic\n"; exit(1);}
		}
		return sol;
	}
	
	template<class U>
	matrix<unsigned int> operator<=(const matrix<U> &mat,U alpha)
	{
		matrix<unsigned int> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			if(mat.data[i] <= alpha) {sol.data[i] = 1;}
			else if(mat.data[i] > alpha) {sol.data[i] = 0;}
			else{std::cout << "Runtime Error\nSomething went wrong with logic\n"; exit(1);}
		}
		return sol;
	}
	
	template<class U>
	matrix<unsigned int> operator<=(U alpha,const matrix<U> &mat)
	{
		matrix<unsigned int> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			if(alpha <= mat.data[i]) {sol.data[i] = 1;}
			else if(alpha > mat.data[i]) {sol.data[i] = 0;}
			else{std::cout << "Runtime Error\nSomething went wrong with logic\n"; exit(1);}
		}
		return sol;
	}
	
	template<class U>
	matrix<unsigned int> operator==(const matrix<U> &mat,U alpha)
	{
		matrix<unsigned int> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			if(alpha == mat.data[i]) {sol.data[i] = 1;}
			else if(alpha != mat.data[i]) {sol.data[i] = 0;}
			else {std::cout << "Runtime Error\nSomething went wrong with the logic\n"; exit(1);}
		}
		return sol;
	}
	
	template<class U>
	matrix<unsigned int> operator==(U alpha,const matrix<U> &mat)
	{
		matrix<unsigned int> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			if(alpha == mat.data[i]) {sol.data[i] = 1;}
			else if(alpha != mat.data[i]) {sol.data[i] = 0;}
			else {std::cout << "Runtime Error\nSomething went wrong with the logic\n"; exit(1);}
		}
		return sol;
	}
	
	template<class U>
	matrix<unsigned int> operator!=(const matrix<U> &mat,U alpha)
	{
		matrix<unsigned int> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			if(alpha == mat.data[i]) {sol.data[i] = 0;}
			else if(alpha != mat.data[i]) {sol.data[i] = 1;}
			else {std::cout << "Runtime Error\nSomething went wrong with the logic\n"; exit(1);}
		}
		return sol;
	}
	
	template<class U>
	matrix<unsigned int> operator!=(U alpha,const matrix<U> &mat)
	{
		matrix<unsigned int> sol(mat.row,mat.col);

		for(int i = 0; i < mat.row*mat.col; i++)
		{
			if(alpha == mat.data[i]) {sol.data[i] = 0;}
			else if(alpha != mat.data[i]) {sol.data[i] = 1;}
			else {std::cout << "Runtime Error\nSomething went wrong with the logic\n"; exit(1);}
		}
		return sol;
	}
	
	template<class T>
	matrix<T> matrix<T>::operator()(const matrix<unsigned int> &mat)
	{
		/*Initally count number of 1s*/
		unsigned int count = 0;

		for(int i = 0; i < row*col; i++)
		{
			if(mat.raw_data()[i] == 0){/*Do nothing*/}
			else if(mat.raw_data()[i] == 1) {count++;}
			else {std::cout << "Runtime Error\nInvalid logical\n"; exit(1);}
		}
		/*Checks the count*/
		if(count == 0) {std::cout << "Runtime Error\nNo element found\n"; exit(1);} else {}
		/*Generate our new matrix*/
		matrix<T> sol;
		/*If statement to preserve the vector(only vectors!) origonal shape*/
		if(mat.rows() == 1) {sol = matrix<T>::shape(1,count);}
		else if(mat.cols() == 1) {sol = matrix<T>::shape(count,1);}
		else {sol = matrix<T>::shape(count,1);}
		//matrix<T> sol(count,1);
		unsigned int index_ = 0;

		for(int i = 0; i < row*col; i++)
		{
			if(mat.raw_data()[i] == 0){/*Do nothing*/}
			else if(mat.raw_data()[i] == 1) {sol.data[index_++] = data[i];}
			else {std::cout << "Runtime Error\nInvalid logical\n"; exit(1);}
		}
		return sol;
	}
	
	/*Find method that returns the indexed values that are true*/
	template<class T>
	matrix<unsigned int> matrix<T>::find(const matrix<unsigned int> &mat)
	{
		/*Counts the size*/
		unsigned int count = 0;

		for(int i = 0; i < row*col; i++)
		{
			if(mat.raw_data()[i] == 0){/*Do nothing*/}
			else if(mat.raw_data()[i] == 1) {count++;}
			else {std::cout << "Runtime Error\nInvalid logical \n"; exit(1);}
		}
		/*Checks the count*/
		if(count == 0) {std::cout << "Runtime Error\nNo element found\n"; exit(1);} else {}
		
		/*Add in indexes*/
		matrix<unsigned int> sol(count,1);
		unsigned int index_ = 0;

		for(int i = 0; i < row*col; i++)
		{
			if(mat.raw_data()[i] == 0){/*Do nothing*/}
			else if(mat.raw_data()[i] == 1) {sol.raw_data()[index_++] = i;}
			else {std::cout << "Runtime Error\nInvalid logical\n"; exit(1);}
		}
		return sol;
	}
	
	template<class T>
	matrix<T> & matrix<T>::set_col(unsigned int c,const matrix<T> &M)
	{
		if(c < 0 || c >= col || M.rows() > this->rows()) {std::cout << "Runtime Error\nInvlaid index\n"; exit(1);} else {}

		for(int i = 0; i < this->rows(); i++)
		{
			data[row*i + c] = M(i);
		}
		return *this; 
	}
	
	template<class T>
	matrix<T> & matrix<T>::set_row(unsigned int r,const matrix<T> &M)
	{
		if(r < 0 || r >= this->cols() || M.cols() > this->cols()) {std::cout << "Runtime Error\nInvlaid index\n"; exit(1);} else {}

		for(int i = 0; i < this->cols(); i++)
		{
			data[r*row + i] = M(i);
		}
		return *this; 
	}
}

/*Name declarations*/
#define mat matrix<double>
#define fmat matrix<float>
#define imat matrix<int>
#define uimat matrix<unsigned int>
#define hpmat matrix<long double>
#define bmat matrix<bool>

/*Complex mats*/
#define cmat matrix<std::complex<double> >
#define cfmat matrix<std::complex<float> >
#define cimat matrix<std::complex<int> >
#define cuimat matrix<std::complex<unsigned int> >
#define chpmat matrix<std::complex<long double> >

/*function defines*/
#define drandu matrix<double>::randu
#define frandu matrix<float>::randu
#define irandu matrix<int>::randu
#define uirandu matrix<unsigned int>::randu
#define hpranu matrix<long double>::randu

#define dzeros matrix<double>::zeros
#define fzeros matrix<float>::zeros
#define izeros matrix<int>::zeros
#define uizeros matrix<unsigned int>::zeros
#define hpzeros matrix<long double>::zeros

#define dones matrix<double>::ones
#define fones matrix<float>::ones
#define iones matrix<int>::ones
#define uiones matrix<unsigned int>::ones
#define hpones matrix<long double>::ones

#define dlinspace matrix<double>::linspace
#define flinspace matrix<float>::linspace
#define ilinspace matrix<int>::linspace
#define uilinspace matrix<unsigned int>::linspace
#define hplinspace matrix<long double>::linspace

#define ddelspace matrix<double>::delspace
#define fdelspace matrix<float>::delspace
#define idelspace matrix<int>::delspace
#define uidelspace matrix<unsigned int>::delspace
#define hpdelspace matrix<long double>::delspace

#define drandn matrix<double>::randn
#define frandn matrix<float>::randn
#define irandn matrix<int>::randn
#define uirandn matrix<unsigned int>::randn
#define hprandn matrix<long double>::randn

#define deye matrix<double>::eye
#define feye matrix<float>::eye
#define ieye matrix<int>::eye
#define uieye matrix<unsigned int>::eye
#define hpeye matrix<long double>::eye

#define dshape matrix<double>::shape
#define fshape matrix<float>::shape
#define ishape matrix<int>::shape
#define uishape matrix<unsigned int>::shape
#define hpshape matrix<long double>::shape

#define ui unsigned int
#define hp long double
#define ll long long




