#ifndef COMPLEX_H
#define COMPLEX_H
#pragma once
class Complex {
private:
	double a, b, r, theta;
public:
	Complex();
	~Complex() {}
	Complex(const Complex& inp);
	Complex& operator=(const Complex& rhs);
	Complex(Complex&& inp);
	Complex& operator=(Complex&& rhs);

	Complex(double A, double B);

	Complex operator+(const Complex& rhs);
	Complex operator-(const Complex& rhs);
	Complex operator*(const Complex& rhs);
	Complex operator/(const Complex& rhs);
	Complex Complex::operator+(double rhs) const;
	Complex Complex::operator-(double rhs) const;
	Complex Complex::operator*(double rhs) const;
	Complex Complex::operator/(double rhs) const;
	friend Complex operator+(double lhs, const Complex& rhs);
	friend Complex operator-(double lhs, const Complex& rhs);
	friend Complex operator*(double lhs, const Complex& rhs);
	friend Complex operator/(double lhs, const Complex& rhs);
	friend Complex conj(const Complex& arg);

	friend Complex exp(const Complex& arg);


	void setPolar(double R, double THTA);
	void setCart(double A, double B);

	double getA() { return a; }
	double getB() { return b; }
	double getR() { return r; }
	double getTheta() { return theta; }
};
#endif