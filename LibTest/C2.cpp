#include "stdud.h"

#include <iostream>
#include "C2.h"
#include "Complex.h"
#include "cR2.h"
#include "VecR2.h"
#include "CUDA_Functions_R2.h"
using namespace std;
C2::C2() {}

C2::~C2() {}

C2::C2(const C2 & inp) {
	a = inp.a; b = inp.b;
	r = inp.r; theta = inp.theta;
}

C2::C2(int M, int N) {
	a = cR2(M, N);
	b = cR2(M, N);
	r = cR2(M, N);
	theta = cR2(M, N);
}

C2::C2(cR2 & A, cR2 & B) {
	//cR2* temp = new cR2();
	if (sameDim(A,B)) {
		//cout << A(1100, 1100) << endl;
		a = A; b = B;
		//a.call2Host();
		//cout << a.getSize() << endl;
		r = sqrt(A*A + B*B);
		theta = atan(B,A);
		//cout << r(1100, 1100) << endl;
	}
}

Complex C2::operator()(int I, int J) {
	return Complex(a(I,J),b(I,J));
}

Complex C2::operator()(double A, double B) {
	return Complex(a(A,B),b(A,B));
}

Complex C2::operator()(VecR2 & pos) {
	return Complex(a(pos.getX(), pos.getY()),b(pos.getX(), pos.getY()));
}

C2 C2::operator+(const C2 & rhs) {
	C2 temp(a + rhs.a, b + rhs.b);
	return temp;
}

C2 C2::operator-(const C2 & rhs) {
	C2 temp(a - rhs.a, b - rhs.b);
	return temp;
}

C2 C2::operator*(const C2 & rhs) {
	C2 temp;
	temp.setPolar(r*rhs.r, theta + rhs.theta);
	return temp;
}

C2 C2::operator/(const C2 & rhs) {
	C2 temp;
	temp.setPolar(r/rhs.r, theta - rhs.theta);
	return temp;
}

C2 C2::operator+(double rhs) const {
	C2 temp(*this);
	temp.a = temp.a + rhs;
	temp.r = sqrt(temp.a*temp.a + temp.b*temp.b);
	temp.theta = atan(temp.b, temp.a);
	return temp;
}

C2 C2::operator-(double rhs) const {
	C2 temp(*this);
	temp.a = temp.a - rhs;
	temp.r = sqrt(temp.a*temp.a + temp.b*temp.b);
	temp.theta = atan(temp.b, temp.a);
	return temp;
}

C2 C2::operator*(double rhs) const {
	C2 temp(*this);
	temp.r = rhs*temp.r;
	temp.a = temp.r*cos(temp.theta);
	temp.b = temp.r*sin(temp.theta);
	return temp;
}

C2 C2::operator/(double rhs) const {
	C2 temp(*this);
	temp.r = temp.r/rhs;
	temp.a = temp.r*cos(temp.theta);
	temp.b = temp.r*sin(temp.theta);
	return temp;
}

void C2::setCart(cR2 & A, cR2 & B) {
	a = A; b = B;
	r = sqrt(a*a + b*b);
	theta = atan(b, a);
}

void C2::setPolar(cR2 & R, cR2 & Theta) {
	r = R; theta = Theta;
	a = r*cos(theta); b = r*sin(theta);
}

void C2::send2Device() {
	a.send2Device();
	b.send2Device();
	r.send2Device();
	theta.send2Device();
}

void C2::call2Host() {
	a.call2Host();
	b.call2Host();
	r.call2Host();
	theta.call2Host();
}

C2 operator+(double lhs, const C2 & rhs) {
	C2 temp(rhs);
	temp.a = lhs + temp.a;
	temp.r = sqrt(temp.a*temp.a + temp.b*temp.b);
	temp.theta = atan(temp.b,temp.a);
	return temp;
}

C2 operator-(double lhs, const C2 & rhs) {
	C2 temp(rhs);
	temp.a = lhs - temp.a;
	temp.r = sqrt(temp.a*temp.a + temp.b*temp.b);
	temp.theta = atan(temp.b, temp.a);
	return temp;
}

C2 operator*(double lhs, const C2 & rhs) {
	C2 temp(rhs);
	temp.r = lhs*rhs.r;
	temp.a = temp.r*cos(temp.theta); temp.b = temp.r*sin(temp.theta);
	return temp;
}

C2 operator/(double lhs, const C2 & rhs) {
	C2 temp(rhs);
	temp.r = lhs/rhs.r;
	temp.theta = (-1.0)*rhs.theta;
	temp.a = temp.r*cos(temp.theta); temp.b = temp.r*sin(temp.theta);
	return temp;
}

C2 operator*(Complex lhs, const cR2 & rhs) {
	C2 temp;
	temp.setCart(lhs.getA()*rhs, lhs.getB()*rhs);
	return temp;
}

C2 conj(const C2 & inp) {
	C2 temp(inp);
	temp.theta = (-1.0)*inp.theta;
	return temp;
}

C2 normalize(const C2 & inp) {
	C2 temp(inp);
	temp.r = normalize(temp.r);
	temp.a = temp.r*cos(temp.theta); temp.b = temp.r*sin(temp.theta);
	return temp;
}
