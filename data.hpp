#pragma once
#include "matrix.h"
#include <fstream>
#include <iostream>
#include <stdlib.h>

namespace nicole
{
	class file_read_write
	{
	private:
		char * file_name;
	public:
		/*Constructors required*/
		file_read_write(void) = default;
		file_read_write(char *);
		~file_read_write(void);
		
		/*Setters*/
		file_read_write & set_name(char *);
		
		/*Read/Write to begining of file*/
		template<typename T> file_read_write & write(const matrix<T>&);
		template<typename T> file_read_write & read(matrix<T>&);
		
		
	};
	
	file_read_write::file_read_write(char * pname) 
	{
		file_name = pname;
	}
	
	file_read_write::~file_read_write(void)
	{
		/*Do nothing*/
	}
	
	file_read_write & file_read_write::set_name(char * pname)
	{
		file_name = pname;
		return *this;
	}	
		
	template<typename T>
	file_read_write & file_read_write::write(const matrix<T> &m)
	{
		/*Check to see if file is open*/
		std::ofstream file;
		file.open(file_name);
		if(!file.is_open()) {std::cerr << "Runtime Error: File not found\n"; exit(1);}
		
		/*Write to the file*/
		for(int i = 0; i < m.rows(); ++i)
		{
			for(int j = 0; j < m.cols(); ++j)
			{
				file << m(i,j) << " ";
			}
			file << "\n";
		}
		file.close();
		return *this;
	}	
	
	template<typename T>
	file_read_write & file_read_write::read(matrix<T> &m) //Make sure m size is decalared beforetime
	{
		/*Check to see if the file exsits*/
		std::ifstream file;
		file.open(file_name);
		if(!file.is_open()) {std::cerr << "Runtime Error: File not found\n"; exit(1);}
		
		/*Reinitalize the matrix with new values*/
		m = matrix<T>::zeros(m.rows(),m.cols());

		/*Read from the file*/
		T temp;
		unsigned int i = 0;
		while(file >> temp)
		{
			m.raw_data()[i++] = temp;
		}

		file.close();
		return *this;
	}	
}