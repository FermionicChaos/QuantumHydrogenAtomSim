#ifndef H_ATOM_V2_H
#define H_ATOM_V2_H
#pragma once
#include <cuda.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "cR2.h"
#include "C2.h"
#include "Complex.h"
#include <SFML/Graphics.hpp>
using namespace std;
using namespace sf;

class H_Atom_v2 {
public:
	VertexArray Im;
	//C2 PreIm;
	Uint8* h_pix;
	Uint8 *d_pix;
	cR2 alpha, phase;

	double tt;
	double *d_psia, *d_psib, *d_phase, *d_R; //Device Superposition wave function
	double *h_psia, *h_psib, *h_phase, *h_R; //Host super position.
	double *h_Ef, *d_Ef; // Eigen function vector storing the multiple eigenfunctions.
	double *h_avec, *h_bvec;
	double *d_avec, *d_bvec; // State vector stored in gpu memory.
	//int length, size; //Length = vector lentgth of the energy eigen states. size = discretized set of real-valued eigenfunctions.
	int I, J, K, h, w, Nh;
	double dx, dy;
	dim3 grid, block;

	H_Atom_v2();
	~H_Atom_v2();

	void StateVec(double* psia, double* psib, cR2 & x, cR2 & y, int N);
	void genData(int RESX, int RESY, cR2& x, cR2& y, int N);
	void genIm(int horizontal, int vertical);

	//void InnerProduct(Vector2i inp);

	void compute(float dt);
	//void update(float dt);
	void draw(RenderWindow& window);
	void ArrayTest();
};



#endif // H_ATOM_V2_H