#ifndef C2_H
#define C2_H
#pragma once
#include "cR2.h"
#include "VecR2.h"
#include "Complex.h"

//Fix me!
class C2 {
private:
	//R2 a, b, r, theta;
	//int m, n;
	//double dx, dy;
	//double x1, x2, y1, y2;
public:
	cR2 a, b, r, theta;
	C2();
	~C2();
	C2(const C2& inp);
	C2(int M, int N);

	C2(cR2& A, cR2& B);

	Complex C2::operator()(int I, int J);
	Complex C2::operator()(double A, double B);
	Complex C2::operator()(VecR2& pos);

	C2 C2::operator+(const C2& rhs);
	C2 C2::operator-(const C2& rhs);
	C2 C2::operator*(const C2& rhs);
	C2 C2::operator/(const C2& rhs);
	
	C2 C2::operator+(double rhs) const;
	C2 C2::operator-(double rhs) const;
	C2 C2::operator*(double rhs) const;
	C2 C2::operator/(double rhs) const;
	friend C2 operator+(double lhs, const C2& rhs);
	friend C2 operator-(double lhs, const C2& rhs);
	friend C2 operator*(double lhs, const C2& rhs);
	friend C2 operator/(double lhs, const C2& rhs);

	friend C2 operator*(Complex lhs, const cR2& rhs);

	friend C2 conj(const C2& inp);
	friend C2 normalize(const C2& inp);

	void setCart(cR2& A, cR2& B);
	void setPolar(cR2& R, cR2& Theta);
	void send2Device();
	void call2Host();

	cR2 getA() { return a; }
	cR2 getB() { return b; }
	cR2 getR() { return r; }
	cR2 getT() { return theta; }
};
#endif // !1
