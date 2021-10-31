//#include "stdafx.h"
#include "stdud.h"
#include "CUDA_Functions_R2.h"
#include "FunctionsR2.h"
using namespace std;

double a0 = 1.0, hbar = 1.0, Z = 1.0;
double e = 1.0, mu = 1.0;


cR2 Plm(int l, int m, cR2 & arg) {
	//cout << l << " " << m << endl;
	arg.call2Host();
	double* ba = arg.getH_ptr();
	int N = arg.getSize();
	cR2 temp(arg);
	temp.call2Host();
	double* bf = temp.getH_ptr();
	for (int i = 0; i < N; ++i) {
		bf[i] = gsl_sf_legendre_Plm(l, m, ba[i]);
	}
	temp.send2Device();
	return temp;
}

cR2 Ylm(int l, int m, cR2 & arg) {
	//cout << l << " " << m << endl;
	arg.call2Host();
	//cout << arg(1100, 1100) << endl;
	double* ba = arg.getH_ptr();
	int N = arg.getSize();
	cR2 temp(arg);
	temp.call2Host();
	//cout << temp(1100, 1100) << endl;
	double* bf = temp.getH_ptr();
	if (m >= 0) {
		for (int i = 0; i < N; ++i) {
			bf[i] = gsl_sf_legendre_sphPlm(l, m, ba[i]);
		}
		temp.send2Device();
		//cout << "Ylm(0.0,0.0) = " << temp(0.0, 0.0) << endl;
		return temp;
	}
	else {
		for (int i = 0; i < N; ++i) {
			bf[i] = pow(-1.0,-m)*gsl_sf_legendre_sphPlm(l, -m, ba[i]);
		}
		temp.send2Device();
		return temp;
	}
}

cR2 Lka(int k, double a, cR2 & arg) {
	//cout << "|" << k << a << ">" << endl;
	arg.call2Host();
	//cout << arg(1100, 1100) << endl;
	double* ba = arg.getH_ptr();
	int N = arg.getSize();
	cR2 temp(arg);
	temp.call2Host();
	//cout << temp(1100, 1100) << endl;
	double* bf = temp.getH_ptr();
	for (int i = 0; i < N; ++i) {
		//L
		bf[i] = gsl_sf_laguerre_n(k, a, ba[i]);
	}
	temp.send2Device();
	return temp;
}

cR2 Rnl(double n, double l, cR2 & arg) {
	cR2 temp;
	cR2 rho = (2 * Z / (n*a0))*arg;
	//rho.call2Host();
	//cout << rho(0.0,0.0) << endl;
	double Cnl = sqrt(4*Z*tgamma(n - l)/(n*n*a0*tgamma(n+l+1)));
	temp = Cnl*pow(rho, l)*Lka(n - l - 1, 2*l + 1, rho)*exp((-0.5)*rho);
	//temp.call2Host();
	//cout << "Rnl(0.0,0.0) = " << temp(0.0, 0.0) << endl;
	return temp;
}

cR2 Psi_nlm(int n, int l, int m, cR2 & xx, cR2 & yy) {
	//Note2self: on a rainy day, design Legendre, Lageurre family for GPU.
	//cout << n << " " << l << " " << m << endl;
	//cR2 argR = sqrt(xx*xx + yy*yy);
	//cR2 argT = yy / sqrt(xx*xx + yy*yy);
	//cR2 temp = Rnl(n, l, argR)*Ylm(l, m, argT);

	cR2 temp = Rnl(n, l, sqrt(xx*xx + yy*yy))*Ylm(l, m, yy / sqrt(xx*xx + yy*yy));
	//temp.call2Host();
	//cout << temp(1100, 1100) << endl;
	return temp;
}


R2 Plm(int l, int m, const R2 & arg) {
	double* ba = arg.getBase();
	int N = arg.getSize();
	R2 temp(arg);
	double* bf = temp.getBase();
	for (int i = 0; i < N; ++i) {
		bf[i] = gsl_sf_legendre_Plm(l, m, ba[i]);
	}
	return temp;
}

R2 Ylm(int l, int m, const R2 & arg) {
	double* ba = arg.getBase();
	int N = arg.getSize();
	R2 temp(arg);
	double* bf = temp.getBase();
	if (m >= 0) {
		for (int i = 0; i < N; ++i) {
			bf[i] = gsl_sf_legendre_sphPlm(l, m, ba[i]);
		}
		return temp;
	}
	else {
		for (int i = 0; i < N; ++i) {
			bf[i] = pow(-1.0, -m)*gsl_sf_legendre_sphPlm(l, -m, ba[i]);
		}
		return temp;
	}
}

R2 Lka(int k, double a, const R2 & arg) {
	double* ba = arg.getBase();
	int N = arg.getSize();
	R2 temp(arg);
	double* bf = temp.getBase();
	for (int i = 0; i < N; ++i) {
		bf[i] = gsl_sf_laguerre_n(k, a, ba[i]);
	}
	return temp;
}

R2 Rnl(double n, double l, const R2 & arg) {
	R2 temp;
	R2 rho = (2 * Z / (n*a0))*arg;
	double Cnl = sqrt(Z*tgamma(n - l) / (n*n*a0*tgamma(n + l + 1)));
	temp = Cnl*pow(rho, l)*Lka(n - l - 1, l, rho)*exp((-1.0)*rho);
	return temp;
}

R2 Psi_nlm(int n, int l, int m, R2 & xx, R2 & yy) {
	R2 temp = Rnl(n, l, sqrt(xx*xx + yy*yy))*Ylm(l, m, yy / sqrt(xx*xx + yy*yy));
	return temp;
}

R2 exp(const R2 & arg) {
	double* ba = arg.getBase();
	int N = arg.getSize();
	R2 temp(arg);
	double* bf = temp.getBase();
	for (int i = 0; i < N; ++i) {
		bf[i] = exp(ba[i]);
	}
	return temp;
}

R2 sqrt(const R2 & arg) {
	double* ba = arg.getBase();
	int N = arg.getSize();
	R2 temp(arg);
	double* bf = temp.getBase();
	for (int i = 0; i < N; ++i) {
		bf[i] = sqrt(ba[i]);
	}
	return temp;
}

R2 pow(const R2 & arg, double exponent) {
	double* ba = arg.getBase();
	int N = arg.getSize();
	R2 temp(arg);
	double* bf = temp.getBase();
	for (int i = 0; i < N; ++i) {
		bf[i] = pow(ba[i], exponent);
	}
	return temp;
}

R2 atan2(const R2 & y, const R2 & x) {
	double* by = y.getBase();
	double* bx = x.getBase();
	int N = y.getSize();
	R2 temp(y);
	double* bf = temp.getBase();
	for (int i = 0; i < N; ++i) {
		bf[i] = atan2(by[i], bx[i]);
	}
	return temp;
}



/*
C2 tevolve(H_Atom& inp) {
	C2 temp;
	float tt = inp.tt;
	double n = static_cast<double>(inp.n);
	double* state_vec = inp.state_vec;
	Complex e;
	int index = 0;
	for (int i = 1; i < inp.n; ++i) {
		for (int j = 0; j < i; ++j) {
			for (int k = -j; k < j; ++k) {
				e = Complex(0.0, 1.0*tt/(2*i*i));
				e = exp(e);
				if (index == 0) {
					temp = e*state_vec[index] * inp.Wave[index];
				}
				temp = temp + e*state_vec[index] * inp.Wave[index];
				index += 1;
			}

		}

	}
	return temp;
}*/

