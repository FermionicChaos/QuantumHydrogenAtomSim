//#include "stdafx.h"
#include "stdud.h"
#include "FunctionsR1.h"

R1 Plm(int l, int m, const R1 & arg) {
	double* ba = arg.getBase();
	int N = arg.getSize();
	R1 temp(arg);
	double* bf = temp.getBase();
	for (int i = 0; i < N; ++i) {
		bf[i] = gsl_sf_legendre_Plm(l, m, ba[i]);
	}
	return temp;
}

R1 Ylm(int l, int m, const R1 & arg) {
	double* ba = arg.getBase();
	int N = arg.getSize();
	R1 temp(arg);
	double* bf = temp.getBase();
	for (int i = 0; i < N; ++i) {
		bf[i] = gsl_sf_legendre_sphPlm(l, m, ba[i]);
	}
	return temp;
}

R1 Lka(int k, double a, const R1 & arg) {
	double* ba = arg.getBase();
	int N = arg.getSize();
	R1 temp(arg);
	double* bf = temp.getBase();
	for (int i = 0; i < N; ++i) {
		bf[i] = gsl_sf_laguerre_n(k, a, ba[i]);
	}
	return temp;
}

