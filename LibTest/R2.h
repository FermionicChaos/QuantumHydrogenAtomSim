#ifndef R2_H
#define R2_H
#pragma once
#include "stdud.h"
#include "VecR2.h"
class R2 {
private:
	double* base;
	int m, n;
	double dx, dy;
	double x1, x2, y1, y2;
public:
	R2(); //Default Constructor
	~R2(); //Destructor
	R2(const R2& inp); //Copy Construct
	R2& operator=(const R2& rhs); //Copy Assignment
	R2(R2&& inp); //Move Construct
	R2& operator=(R2&& rhs); //Move Assignment
	R2(int M, int N); //Reserve Memory

	/*Specialized Initializers*/
	R2(int k, double X1, double X2, int M, double Y1, double Y2, int N); //Create set domain.
	friend std::ostream &operator<<(std::ostream &os, const R2& rhs);
	
	/*Access*/
	double R2::operator()(int I, int J);
	double R2::operator()(double X, double Y);
	double R2::operator()(VecR2& pos);
	//double R2::operator()(const Vec2& V);

	R2 operator+(const R2& g);
	R2 operator-(const R2& g);
	R2 operator*(const R2& g);
	R2 operator/(const R2& g);

	R2 R2::operator+(double rhs);
	R2 R2::operator-(double rhs);
	R2 R2::operator*(double rhs);
	R2 R2::operator/(double rhs);
	friend R2 operator+(double lhs, const R2& rhs);
	friend R2 operator-(double lhs, const R2& rhs);
	friend R2 operator*(double lhs, const R2& rhs);
	friend R2 operator/(double lhs, const R2& rhs);

	friend R2 sin(const R2& arg);
	friend R2 cos(const R2& arg);
	friend R2 tan(const R2& arg);
	friend R2 asin(const R2& arg);
	friend R2 acos(const R2& arg);
	friend R2 atan(const R2& arg);
	//friend R2 atan(const R2& argy, const R2& argx);
	
	friend R2 sinh(const R2& arg);
	friend R2 cosh(const R2& arg);
	friend R2 tanh(const R2& arg);

	friend R2 d_dx(const R2& f, int Ver);
	friend R2 d_dy(const R2& f, int Ver);
	friend R2 normalize(R2& f);

	int getIndex1() const { return m; }
	int getIndex2() const { return n; }
	double getDX() const { return dx; }
	double getDY() const { return dy; }
	double getX1() const { return x1; }
	double getX2() const { return x2; }
	double getY1() const { return y1; }
	double getY2() const { return y2; }
	double* getBase() const { return base; }
	int getMEM() const { return sizeof(double)*m*n; }
	int getSize() const { return m*n; }
	double getMax() const;

	void setIndex1(int M) { m = M; }
	void setIndex2(int N) { n = N; }
	void setDX(double DX) { dx = DX; }
	void setDY(double DY) { dy = DY; }
	void setX1(double X1) { x1 = X1; }
	void setX2(double X2) { x2 = X2; }
	void setY1(double Y1) { y1 = Y1; }
	void setY2(double Y2) { y2 = Y2; }
	void setBase(double* B) { base = B; }
	bool sameDim(const R2& f, const R2& g);
};
#endif