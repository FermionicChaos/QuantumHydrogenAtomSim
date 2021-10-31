#ifndef CUDAR2_H
#define CUDAR2_H
#pragma once
#include <cuda.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "stdud.h"
#include "VecR2.h"

// Thread number set at 512 per block.

class cudaR2 {
private:
	double* ptr;
	int m, n;
	double dx, dy;
	double x1, x2, y1, y2;
	dim3 grid, block;
public:
	cudaR2(); //Default Constructor
	~cudaR2(); //Destructor
	cudaR2(const cudaR2& inp); //Copy Construct
	cudaR2& operator=(const cudaR2& rhs); //Copy Assignment
	cudaR2(cudaR2&& inp); //Move Construct
	cudaR2& operator=(cudaR2&& rhs); //Move Assignment
	cudaR2(int M, int N); //Reserve Memory

	cudaR2(int k, double X1, double X2, int M, double Y1, double Y2, int N); //Create set domain.

	double cudaR2::operator()(int I, int J);
	double cudaR2::operator()(double X, double Y);
	double cudaR2::operator()(VecR2& pos);

	cudaR2 operator+(const cudaR2& rhs);
	cudaR2 operator-(const cudaR2& rhs);
	cudaR2 operator*(const cudaR2& rhs);
	cudaR2 operator/(const cudaR2& rhs);
	
	cudaR2 cudaR2::operator+(double rhs);
	cudaR2 cudaR2::operator-(double rhs);
	cudaR2 cudaR2::operator*(double rhs);
	cudaR2 cudaR2::operator/(double rhs);

	friend cudaR2 operator+(double lhs, const cudaR2& rhs);
	friend cudaR2 operator-(double lhs, const cudaR2& rhs);
	friend cudaR2 operator*(double lhs, const cudaR2& rhs);
	friend cudaR2 operator/(double lhs, const cudaR2& rhs);

	int getIndex1() const { return m; }
	int getIndex2() const { return n; }
	double getDX() const { return dx; }
	double getDY() const { return dy; }
	double getX1() const { return x1; }
	double getX2() const { return x2; }
	double getY1() const { return y1; }
	double getY2() const { return y2; }
	double* getBase() const { return ptr; }
	int getMEM() const { return sizeof(double)*m*n; }
	int getSize() const { return m*n; }
	dim3 getG() const { return grid; }
	dim3 getB() const { return block; }
	//double getMax() const;

	void setIndex1(int M) { m = M; }
	void setIndex2(int N) { n = N; }
	void setDX(double DX) { dx = DX; }
	void setDY(double DY) { dy = DY; }
	void setX1(double X1) { x1 = X1; }
	void setX2(double X2) { x2 = X2; }
	void setY1(double Y1) { y1 = Y1; }
	void setY2(double Y2) { y2 = Y2; }
	void setBase(double* B) { ptr = B; }
	void setG(dim3 G) { grid = G; }
	void setB(dim3 B) { block = B; }
	friend bool sameDim(const cudaR2& f, const cudaR2& g);
};
#endif // CUDAR2_H