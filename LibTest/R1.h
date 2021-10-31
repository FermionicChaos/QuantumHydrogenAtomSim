#ifndef R1_H
#define R1_H
#pragma once
//Describe functions along discrete lines.
//Include size, ptr.
//Add Scalar operations.
#include <iostream>
class R1 {
private:
	double* base; //Base address on heap.
	int n; //Elements of the
	double dx; //interval
	double x1, x2;
public:

	//Constructors x2 > x1
	R1(); //Default Constructor (No Reserve)
	~R1(); //Destructor
	R1(const R1& inp); // Copy Constructor R1 yy = xx;
	R1& operator=(const R1& rhs); // yy = xx
	R1(R1&& inp);
	R1& operator=(R1&& rhs);
	R1(int size); //Reserve space
	
	R1(double a, double b, int size); // Extrapolate dx
	
	friend std::ostream &operator<<(std::ostream &os, const R1& rhs);
	double R1::operator()(double x0);
	double R1::operator()(int index);
	//Operators
	R1 operator+(const R1& g);
	R1 operator-(const R1& g);
	R1 operator*(const R1& g);
	R1 operator/(const R1& g);
	R1 R1::operator+(double g);
	R1 R1::operator-(double g);
	R1 R1::operator*(double g);
	R1 R1::operator/(double g);
	friend R1 operator+(double lhs, const R1& rhs);
	friend R1 operator-(double lhs, const R1& rhs);
	friend R1 operator*(double lhs, const R1& rhs);
	friend R1 operator/(double lhs, const R1& rhs);


	//Member Functions
	friend R1 cos(const R1& x);
	friend R1 sin(const R1& x);
	friend R1 tan(const R1& x);
	friend R1 acos(const R1& x);
	friend R1 asin(const R1& x);
	friend R1 atan(const R1& x);
	friend R1 atan2(const R1& y, const R1& x);

	friend R1 sinh(const R1& x);
	friend R1 cosh(const R1& x);
	friend R1 tanh(const R1& x);
	friend R1 asinh(const R1& x);
	friend R1 acosh(const R1& x);
	friend R1 atanh(const R1& x);
	friend R1 exp(const R1& x);
	friend R1 ln(const R1& x);
	friend R1 pow(double base, const R1& exp);
	friend R1 pow(const R1& base, double exp);
	friend R1 pow(const R1& base, const R1& exp);


	friend R1 derivative(const R1& x, int Ver);
	friend double integralD(double a, double b, const R1& f);
	friend double int_dx(int i1, int i2, const R1& f);



	//Access and Set
	int getIndex(double p);
	double* getBase() const { return base; }
	int getSize() const { return n; }
	double getDX() const { return dx; }
	double getX1() const { return x1; }
	double getX2() const { return x2; }

	void setBase(double* ptr) { base = ptr; }
	void setSize(int N) { n = N; }
	void setDX(double DX) { dx = DX; }
	void setX1(double X1) { x1 = X1; }
	void setX2(double X2) { x2 = X2; }

};
#endif