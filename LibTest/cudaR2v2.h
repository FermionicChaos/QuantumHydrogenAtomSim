#ifndef CUDAR2V2_H
#define CUDAR2V2_H
#include <cuda.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "stdud.h"
#include "VecR2.h"
#pragma once

class cudaR2v2 {
private:
	double* h_ptr;
	double *d_ptr;
	int m, n;
	double dx, dy; // Length element 0.1 0.6 0.3 -0.1
	double x1, x2, y1, y2; //Bounds of the number plane.
	dim3 grid, block;
public:
	cudaR2v2(); //Default Constructor
	~cudaR2v2(); //Destructor
	cudaR2v2(const cudaR2v2& inp); //Copy Construct
	cudaR2v2& operator=(const cudaR2v2& rhs); //Copy Assignment
	cudaR2v2(cudaR2v2&& inp); //Move Construct
	cudaR2v2& operator=(cudaR2v2&& rhs); //Move Assignment
	cudaR2v2(int M, int N); //Reserve Memory

	cudaR2v2(int k, double X1, double X2, int M, double Y1, double Y2, int N); //Create & set domain.

	double cudaR2v2::operator()(int I, int J);
	double cudaR2v2::operator()(double X, double Y);
	double cudaR2v2::operator()(VecR2& pos);
	//Fully operational + - * /
	cudaR2v2 operator+(const cudaR2v2& rhs);
	cudaR2v2 operator-(const cudaR2v2& rhs);
	cudaR2v2 operator*(const cudaR2v2& rhs);
	cudaR2v2 operator/(const cudaR2v2& rhs);
	//Yeah, I know, poor notation. At least there is easy to read symmettry in the
	//notation.
	friend cudaR2v2 operator+(const cudaR2v2& lhs, double rhs);
	friend cudaR2v2 operator-(const cudaR2v2& lhs, double rhs);
	friend cudaR2v2 operator*(const cudaR2v2& lhs, double rhs);
	friend cudaR2v2 operator/(const cudaR2v2& lhs, double rhs);
	//Fuck these operations below!
	friend cudaR2v2 operator+(double lhs, const cudaR2v2& rhs);
	friend cudaR2v2 operator-(double lhs, const cudaR2v2& rhs);
	friend cudaR2v2 operator*(double lhs, const cudaR2v2& rhs);
	friend cudaR2v2 operator/(double lhs, const cudaR2v2& rhs);

	//These operators have NOT been configured to the GPU.
	friend cudaR2v2 d_dx(const cudaR2v2& f, int Ver);
	friend cudaR2v2 d_dy(const cudaR2v2& f, int Ver);
	friend cudaR2v2 normalize(cudaR2v2& f);
	
	int getIndex1() const { return m; }
	int getIndex2() const { return n; }
	double getDX() const { return dx; }
	double getDY() const { return dy; }
	double getX1() const { return x1; }
	double getX2() const { return x2; }
	double getY1() const { return y1; }
	double getY2() const { return y2; }
	double* getH_ptr() const { return h_ptr; }
	double* getD_ptr() const { return d_ptr; }
	int getMEM() const { return sizeof(double)*m*n; }
	int getSize() const { return m*n; }
	double getMax() const;
	dim3 getG() const { return grid; }
	dim3 getB() const { return block; }

	void setIndex1(int M) { m = M; }
	void setIndex2(int N) { n = N; }
	void setDX(double DX) { dx = DX; }
	void setDY(double DY) { dy = DY; }
	void setX1(double X1) { x1 = X1; }
	void setX2(double X2) { x2 = X2; }
	void setY1(double Y1) { y1 = Y1; }
	void setY2(double Y2) { y2 = Y2; }
	void setH_ptr(double* B) { h_ptr = B; }
	void setD_ptr(double* B) { d_ptr = B; }
	void setG(dim3 G) { grid = G; }
	void setB(dim3 B) { block = B; }

	void send2Device();
	void call2Host();
	friend bool sameDim(const cudaR2v2& f, const cudaR2v2& g);
};


#endif // CUDAR2V2_H
