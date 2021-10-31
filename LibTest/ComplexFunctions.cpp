//#include "stdafx.h"
#include "stdud.h"
#include "ComplexFunctions.h"
#include "Complex.h"

Complex sin(Complex & arg) {
	Complex temp;
	double a0 = arg.getA();
	double b0 = arg.getB();
	double a, b;
	a = sin(a0)*cosh(b0);
	b = sinh(b0)*cos(a);
	temp.setCart(a, b);
	return temp;
}

Complex cos(Complex & arg) {
	Complex temp;
	double a0 = arg.getA();
	double b0 = arg.getB();
	double a, b;
	a = cos(a0)*cosh(b0);
	b = sinh(b0)*sin(a);
	temp.setCart(a, b);
	return temp;
}

Complex tan(Complex & arg) {
	Complex temp;
	temp = sin(arg) / cos(arg);
	return temp;
}
