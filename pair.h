#pragma once
#include <iostream>
#include <stdlib.h>

namespace nicole	
{

template<class T,class U>
class pair
{
private:
	T* lhs;
	U* rhs;
public:
	/*Constructors and destructors*/
	pair(void);
	pair(const T&,const U&);
	pair(const pair<T,U>&);
	~pair(void);
	/*Acessors to left and right*/
	T & left(void);
	U & right(void);
	/*MAthmatical operators*/
	pair<T,U> & operator=(const pair<T,U>&);
	/*Overload to print out results*/
	template<class R,class Y> friend std::ostream & operator<<(std::ostream&,const pair<R,Y>&);
	/*Some functionsal functions*/
	static pair<U,T> swap(const pair<T,U>&);
};

template<class T,class U>
pair<T,U>::pair(void)
{
	lhs = new T;
	rhs = new U;
	if(!lhs || !rhs) {std::cout << "Runtime Error\nCannot allocate memeory\n"; exit(1);} else {}
}

template<class T,class U>
pair<T,U>::pair(const T &l,const U &r)
{
	lhs = new T;
	rhs = new U;
	if(!lhs || !rhs) {std::cout << "Runtime Error\nCannot allocate memeory\n"; exit(1);} else {}
	
	*lhs = l;
	*rhs = r;
}
template<class T,class U>
pair<T,U>::pair(const pair<T,U> &p)
{
	lhs = new T;
	rhs = new U;
	if(!lhs || !rhs) {std::cout << "Runtime Error\nCannot allocate memeory\n"; exit(1);} else {}
	
	*lhs = *(p.lhs);
	*rhs = *(p.rhs);
}

template<class T,class U>
pair<T,U>::~pair(void)
{
	delete rhs;
	delete lhs;
}

template<class T,class U>
T & pair<T,U>::left(void)
{
	return *lhs;
}

template<class T,class U>
U & pair<T,U>::right(void)
{
	return *rhs;
}

template<class T,class U>
pair<T,U> & pair<T,U>::operator=(const pair<T,U> &p)
{
	*rhs = *(p.rhs);
	*lhs = *(p.lhs);
	return *this;
}

template<class R,class Y>
std::ostream & operator<<(std::ostream &s,const pair<R,Y> &p)
{
	s << "(" << (*(p.lhs)) << "," << (*(p.rhs)) << ")";
	return s;
}

template<class T,class U>
pair<U,T> pair<T,U>::swap(const pair<T,U> &p)
{
	pair<U,T> swapped;
	swapped.left() = *(p.rhs);
	swapped.right() = *(p.lhs);
	return swapped;
}








}