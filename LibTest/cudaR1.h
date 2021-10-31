#ifndef CUDAR1_H
#define CUDAR1_H
#pragma once
#include <cuda.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
class cudaR1 {
private:
	double* ptr;
	int m;
	double x1, x2, dx;
	int grid, block;
public:
	cudaR1(); //Default Constructor (No Reserve)
	~cudaR1(); //Destructor
	cudaR1(const cudaR1& inp); // Copy Constructor cudaR1 yy = xx;
	cudaR1& operator=(const cudaR1& rhs); // yy = xx
	cudaR1(cudaR1&& inp);
	cudaR1& operator=(cudaR1&& rhs);
	cudaR1(int size); //Reserve space

	cudaR1(double a, double b, int res);

	//General Access
	double cudaR1::operator()(double p0);
	double cudaR1::operator()(int i);

	//Fully operational Arithmatic operators
	cudaR1 operator+(const cudaR1& rhs);
	cudaR1 operator-(const cudaR1& rhs);
	cudaR1 operator*(const cudaR1& rhs);
	cudaR1 operator/(const cudaR1& rhs);

	cudaR1 cudaR1::operator+(double rhs);
	cudaR1 cudaR1::operator-(double rhs);
	cudaR1 cudaR1::operator*(double rhs);
	cudaR1 cudaR1::operator/(double rhs);

	friend cudaR1 operator+(double lhs, const cudaR1& rhs);
	friend cudaR1 operator-(double lhs, const cudaR1& rhs);
	friend cudaR1 operator*(double lhs, const cudaR1& rhs);
	friend cudaR1 operator/(double lhs, const cudaR1& rhs);

	//friend cudaR1 d_dx(const cudaR1& f, int ver);

	double* getBase() const { return ptr; }
	int getSize() const { return m; }
	double getDX() const { return dx; }
	double getX1() const { return x1; }
	double getX2() const { return x2; }
	int getG() const { return grid; }
	int getB() const { return block; }

	void setBase(double* Np) { ptr = Np; }
	void setSize(int M) { m = M; }
	void setDX(double DX) { dx = DX; }
	void setX1(double X1) { x1 = X1; }
	void setX2(double X2) { x2 = X2; }
	void InitGB();
};


#endif // CUDAR1_H