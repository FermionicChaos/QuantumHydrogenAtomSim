//#include "stdafx.h"
#include "stdud.h"
#include "Complex.h"

const double pi = 3.141592653589793238;
const double e = 2.718281828459045235;

Complex::Complex() {
	a = 0.0; b = 0.0; r = 0.0, theta = 0.0;
}

Complex::Complex(const Complex & inp) {
	this->a = inp.a; this->b = inp.b;
	this->r = inp.r; this->theta = inp.theta;
}

Complex & Complex::operator=(const Complex & rhs) {
	this->a = rhs.a; this->b = rhs.b;
	this->r = rhs.r; this->theta = rhs.theta;
	return *this;
}

Complex::Complex(Complex && inp) {
	this->a = inp.a; this->b = inp.b;
	this->r = inp.r; this->theta = inp.theta;
}

Complex & Complex::operator=(Complex && rhs) {
	this->a = rhs.a; this->b = rhs.b;
	this->r = rhs.r; this->theta = rhs.theta;
	return *this;
}

Complex::Complex(double A, double B) {
	a = A; b = B;
	r = sqrt(A*A + B*B);
	theta = atan2(B, A);
}

Complex Complex::operator+(const Complex & rhs) {
	Complex temp(a + rhs.a, b + rhs.b);
	return temp;
}

Complex Complex::operator-(const Complex & rhs) {
	Complex temp(a - rhs.a, b - rhs.b);
	return temp;
}

Complex Complex::operator*(const Complex & rhs) {
	Complex temp;
	temp.setPolar(r*rhs.r, theta + rhs.theta);
	return temp;
}

Complex Complex::operator/(const Complex & rhs) {
	Complex temp;
	if (rhs.r != 0.0) {
		temp.setPolar(r / rhs.r, theta - rhs.theta);
		return temp;
	}
	else {
		return temp;
	}
}

Complex Complex::operator+(double rhs) const {
	Complex temp(a + rhs, b);
	return temp;
}

Complex Complex::operator-(double rhs) const {
	Complex temp(a - rhs, b);
	return temp;
}

Complex Complex::operator*(double rhs) const {
	Complex temp;
	temp.setPolar(r*rhs,theta);
	return temp;
}

Complex Complex::operator/(double rhs) const {
	Complex temp;
	if (rhs != 0) {
		temp.setPolar(r / rhs, theta);
		return temp;
	}
	else {
		temp.setPolar(1.0, 0.0);
		return temp;
	}
}

void Complex::setPolar(double R, double THTA) {
	if (R == abs(R)) {
		r = R; theta = THTA;
		a = r*cos(theta); b = r*sin(theta);
	}
	else {
		r = -R; theta = THTA + pi;
		a = r*cos(theta); b = r*sin(theta);
	}
}

void Complex::setCart(double A, double B) {
	a = A; b = B;
	r = sqrt(A*A + B*B);
	theta = atan2(B, A);
}

Complex operator+(double lhs, const Complex & rhs) {
	Complex temp;
	temp.setCart(lhs + rhs.a, rhs.b);
	return temp;
}

Complex operator-(double lhs, const Complex & rhs) {
	Complex temp;
	temp.setCart(lhs - rhs.a, -rhs.b);
	return temp;
}

Complex operator*(double lhs, const Complex & rhs) {
	Complex temp;
	temp.setPolar(lhs*rhs.r, rhs.theta);
	return temp;
}

Complex operator/(double lhs, const Complex & rhs) {
	Complex temp;
	temp.setPolar(lhs/rhs.r, -rhs.theta);
	return temp;
}

Complex conj(const Complex & arg) {
	Complex temp;
	temp.setPolar(arg.r, -arg.theta);
	return temp;
}

Complex exp(const Complex & arg) {
	Complex g;
	g.setPolar(exp(arg.a),arg.b);
	return g;
}
