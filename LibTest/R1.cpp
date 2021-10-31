//#include "stdafx.h"
#include "stdud.h"
#include "R1.h"
using namespace std;


/*DO NOT ALLOW TWO OBEJECTS TO SHARE THE SAME POINTER
ADDRESS! CREATE A NEW POINTER AND ALLOCATE MEMORY FOR EACH ASSIGNMENT.

Designed to created numerical sets and used in other functions as well.

Cases Tested:
	-Assignment:
		y = x;
	-Copy:
		R1 y(x)
		R1 y = x;
	-Destructor:
		Operational, R1 destructs when exiting local region.
		Operational, R1* destructs on commands, not on scope.
	-Allocation:
		Dynamic allocation is ready.
	-Operator(+,-,*,/)
		Works R1 z = (x + y);
		Works, rvalues accounted for z = (x + y);

*/

R1::R1() {
	//Default Constructor
	base = nullptr;
	n = 0; dx = 0; x1 = 0; x2 = 0;
	//cout << "Null Constructor: ADDRESS -> " << base << endl;
}

R1::~R1() {
	//Destructor
	//cout << "MEMORY De-Allocated: ADDRESS -> " << base << endl;
	delete[] base;
}

R1::R1(const R1& rhs) {
	//Copy Constructor SEPERATE POINTER VALUES FOR DIFFERENT OBJECTS
	//cout << "Copy Constructor" << endl;
	n = rhs.getSize(); dx = rhs.getDX();
	x1 = rhs.getX1(); x2 = rhs.getX2();
	double* BR = rhs.base;
	base = new double[n];
	//cout << "COPIED From: ADDRESS -> " << rhs.getBase() << endl;
	//cout << "COPIED To: -> " << base << endl;
	for (int i = 0; i < n; ++i) {
		base[i] = BR[i];
	}
	//cout << "END: COPY Constructor Complete" << endl;
}

R1& R1::operator=(const R1& rhs) {
	//cout << "Operation: ASSIGNMENT" << endl;
	//cout << "LEFT ADDRESS -> " << base << endl;
	//cout << "RIGHT ADDRESS -> " << rhs.getBase() << endl;
	//cout << rhs << endl;
	if ((base != rhs.getBase())) {
		this->n = rhs.getSize(); this->dx = rhs.getDX();
		this->x1 = rhs.getX1(); this->x2 = rhs.getX2();
		delete[] this->base;
		base = new double[n];
		for (int i = 0; i < n; ++i) {
			base[i] = rhs.getBase()[i];
		}
		this->base = base;
		return *this;
	}
	else {
		return *this;
	}
}

R1::R1(R1&& inp) {
	base = inp.base;
	inp.base = nullptr;
	n = inp.n; dx = inp.dx;
	x1 = inp.x1; x2 = inp.x2;
}

R1& R1::operator=(R1&& rhs) {
	delete[] base;
	base = rhs.base;
	n = rhs.n; dx = rhs.dx;
	x1 = rhs.x1; x2 = rhs.x2;
	rhs.base = nullptr;
	return *this;
}

R1::R1(int size) {
	//object shell allocation.
	n = size; x1 = 0; x2 = 0; dx = 0;
	base = new double[n];
	//cout << "MEMORY Allocated: ADDRESS -> " << base << endl;
}

R1::R1(double a, double b, int size) {
	n = size; x1 = a; x2 = b;
	double N = static_cast<double>(size);
	dx = (b - a) / N;
	double j = 0;
	double* temp = new double[n];
	for (int i = 0; i < n; ++i) {
		temp[i] = a + j*dx;
		j = j + 1;
	}
	base = temp;
}

std::ostream &operator<<(std::ostream &os, R1 const &rhs) {
	return os << "DOMAIN: " << " [" << rhs.getX1() << "," 
		<< rhs.getX2() << "] Cardinality: " << rhs.getSize() << " dx = " << rhs.getDX();
}

double R1::operator()(double p) {
	//Works on Domain
	if ((p >= x1) && (p <= x2)) {
		double N = (p - x1) / dx;
		int N2 = static_cast<int>(N);
		return base[N2];
	}
	else {
		cout << "NOT AN ELEMENT!" << endl;
		return 0.0;
	}
}

double R1::operator()(int index) {
	if((index < n)||(index > 0)) {
		return base[index];
	}
	else {
		cout << "NOT AN ELEMENT" << endl;
		return 0.0;
	}
}

R1 R1::operator+(const R1& g) {
	//cout << "Operation: ADDITION" << endl;
	if( (n == g.getSize()) && (x1 = g.getX1()) && (x2 = g.getX2()) ) {
		R1 temp(n);
		double* BOUT = temp.getBase();
		double* BL = base;
		double* BR = g.getBase();
		for (int i = 0; i < n; ++i) {
			BOUT[i] = BL[i] + BR[i];
		}
		temp.setDX(dx);
		temp.setX1(x1); temp.setX2(x2);
		//cout << "ADDITION Success: ADDRESS -> " << BOUT << endl; 
		return temp;
	}
	else {
		cout << "INCOMPATIBLE SETS!" << endl;
		R1 temp;
		return temp;
	}
}

R1 R1::operator-(const R1& g) {
	//cout << "Operation: SUBTRACTION" << endl;
	if ((n == g.getSize()) && (x1 = g.getX1()) && (x2 = g.getX2())) {
		R1 temp(n);
		double* BOUT = temp.getBase();
		double* BL = base;
		double* BR = g.getBase();
		for (int i = 0; i < n; ++i) {
			BOUT[i] = BL[i] - BR[i];
		}
		temp.setDX(dx);
		temp.setX1(x1); temp.setX2(x2);
		return temp;
	}
	else {
		cout << "INCOMPATIBLE SETS!" << endl;
		R1 temp;
		return temp;
	}
}

R1 R1::operator*(const R1& g) {
	//cout << "Operation: SUBTRACTION" << endl;
	if ((n == g.getSize()) && (x1 = g.getX1()) && (x2 = g.getX2())) {
		R1 temp(n);
		double* BOUT = temp.getBase();
		double* BL = base;
		double* BR = g.getBase();
		for (int i = 0; i < n; ++i) {
			BOUT[i] = BL[i] * BR[i];
		}
		temp.setDX(dx);
		temp.setX1(x1); temp.setX2(x2);
		return temp;
	}
	else {
		cout << "INCOMPATIBLE SETS!" << endl;
		R1 temp;
		return temp;
	}
}

R1 R1::operator/(const R1& g) {
	//cout << "Operation: SUBTRACTION" << endl;
	if ((n == g.getSize()) && (x1 = g.getX1()) && (x2 = g.getX2())) {
		R1 temp(n);
		double* BOUT = temp.getBase();
		double* BL = base;
		double* BR = g.getBase();
		for (int i = 0; i < n; ++i) {
			if (BR[i] != 0) {
				BOUT[i] = BL[i] / BR[i];
			}
			else {
				BOUT[i] = NAN;
			}
		}
		temp.setDX(dx);
		temp.setX1(x1); temp.setX2(x2);
		return temp;
	}
	else {
		cout << "INCOMPATIBLE SETS!" << endl;
		R1 temp;
		return temp;
	}
}

R1 R1::operator+(double g) {
		R1 temp(n);
		double* BOUT = temp.getBase();
		double* BL = base;
		for (int i = 0; i < n; ++i) {
			BOUT[i] = BL[i] + g;
		}
		temp.setDX(dx);
		temp.setX1(x1); temp.setX2(x2);
		return temp;
}

R1 R1::operator-(double g) {
	R1 temp(n);
	double* BOUT = temp.getBase();
	double* BL = base;
	for (int i = 0; i < n; ++i) {
		BOUT[i] = BL[i] - g;
	}
	temp.setDX(dx);
	temp.setX1(x1); temp.setX2(x2);
	return temp;
}

R1 R1::operator*(double g) {
	R1 temp(n);
	double* BOUT = temp.getBase();
	double* BL = base;
	for (int i = 0; i < n; ++i) {
		BOUT[i] = BL[i] * g;
	}
	temp.setDX(dx);
	temp.setX1(x1); temp.setX2(x2);
	return temp;
}

R1 R1::operator/(double g) {
	R1 temp(n);
	double* BOUT = temp.getBase();
	double* BL = base;
	for (int i = 0; i < n; ++i) {
		if (BOUT[i] != 0) {
			BOUT[i] = BL[i] / g;
		}
		else {
			BOUT[i] = NAN;
		}
	}
	temp.setDX(dx);
	temp.setX1(x1); temp.setX2(x2);
	return temp;
}


R1 operator+(double lhs, const R1 & rhs) {
	R1 temp(rhs.n);
	double* BOUT = temp.getBase();
	double* BR = rhs.base;
	int N = rhs.n;
	for (int i = 0; i < N; ++i) {
		BOUT[i] = lhs + BR[i];
	}
	temp.setDX(rhs.dx);
	temp.setX1(rhs.x1); temp.setX2(rhs.x2);
	return temp;
}

R1 operator-(double lhs, const R1 & rhs) {
	R1 temp(rhs.n);
	double* BOUT = temp.getBase();
	double* BR = rhs.base;
	int N = rhs.n;
	for (int i = 0; i < N; ++i) {
		BOUT[i] = lhs - BR[i];
	}
	temp.setDX(rhs.dx);
	temp.setX1(rhs.x1); temp.setX2(rhs.x2);
	return temp;
}

R1 operator*(double lhs, const R1 & rhs) {
	R1 temp(rhs.n);
	double* BOUT = temp.getBase();
	double* BR = rhs.base;
	int N = rhs.n;
	for (int i = 0; i < N; ++i) {
		BOUT[i] = lhs * BR[i];
	}
	temp.setDX(rhs.dx);
	temp.setX1(rhs.x1); temp.setX2(rhs.x2);
	return temp;
}

R1 operator/(double lhs, const R1 & rhs) {
	R1 temp(rhs.n);
	double* BOUT = temp.getBase();
	double* BR = rhs.base;
	int N = rhs.n;
	for (int i = 0; i < N; ++i) {
		if (BR[i] != 0) {
			BOUT[i] = lhs / BR[i];
		}
		else {
			BOUT[i] = NAN;
		}
	}
	temp.setDX(rhs.dx);
	temp.setX1(rhs.x1); temp.setX2(rhs.x2);
	return temp;
}

R1 cos(const R1& x) {
	R1 temp(x);
	int N = x.getSize();
	for (int i = 0; i < N; ++i) {
		temp.base[i] = cos( x.base[i] );
	}
	return temp;
}

R1 sin(const R1 & x) {
	R1 Temp(x);
	int N = x.getSize();
	for (int i = 0; i < N; ++i) {
		Temp.base[i] = sin(x.base[i]);
	}
	return Temp;
}

R1 tan(const R1 & x) {
	R1 Temp(x);
	int N = x.getSize();
	for (int i = 0; i < N; ++i) {
		Temp.base[i] = tan(x.base[i]);
	}
	return Temp;
}

R1 acos(const R1 & x) {
	R1 Temp(x);
	int N = x.getSize();
	for (int i = 0; i < N; ++i) {
		Temp.base[i] = acos(x.base[i]);
	}
	return Temp;
}

R1 asin(const R1 & x) {
	R1 Temp(x);
	int N = x.getSize();
	for (int i = 0; i < N; ++i) {
		Temp.base[i] = asin(x.base[i]);
	}
	return Temp;
}

R1 atan(const R1 & x) {
	R1 Temp(x);
	int N = x.getSize();
	for (int i = 0; i < N; ++i) {
		Temp.base[i] = atan(x.base[i]);
	}
	return Temp;
}

R1 atan2(const R1 & y, const R1 & x) {
	R1 Temp(x);
	if (x.getSize() == y.getSize()) {
		int N = x.getSize();
		for (int i = 0; i < N; ++i) {
			Temp.base[i] = atan2(y.base[i], x.base[i]);
		}
		cout << "NOT THE SAME DIMENSIONS!" << endl;
		return Temp;
	}
	else {
		return Temp;
	}
}

R1 sinh(const R1 & x) {
	R1 Temp(x);
	int N = x.getSize();
	for (int i = 0; i < N; ++i) {
		Temp.base[i] = sinh(x.base[i]);
	}
	return Temp;
}

R1 cosh(const R1 & x) {
	R1 Temp(x);
	int N = x.getSize();
	for (int i = 0; i < N; ++i) {
		Temp.base[i] = cosh(x.base[i]);
	}
	return Temp;
}

R1 tanh(const R1 & x) {
	R1 Temp(x);
	int N = x.getSize();
	for (int i = 0; i < N; ++i) {
		Temp.base[i] = tanh(x.base[i]);
	}
	return Temp;
}

R1 asinh(const R1 & x) {
	R1 Temp(x);
	int N = x.getSize();
	for (int i = 0; i < N; ++i) {
		Temp.base[i] = asinh(x.base[i]);
	}
	return Temp;
}

R1 acosh(const R1 & x) {
	R1 Temp(x);
	int N = x.getSize();
	for (int i = 0; i < N; ++i) {
		Temp.base[i] = acosh(x.base[i]);
	}
	return Temp;
}

R1 atanh(const R1 & x) {
	R1 Temp(x);
	int N = x.getSize();
	for (int i = 0; i < N; ++i) {
		Temp.base[i] = atanh(x.base[i]);
	}
	return Temp;
}

R1 exp(const R1 & x) {
	R1 Temp(x);
	int N = x.getSize();
	for (int i = 0; i < N; ++i) {
		Temp.base[i] = exp(x.base[i]);
	}
	return Temp;
}

R1 ln(const R1 & x) {
	R1 Temp(x);
	int N = x.getSize();
	for (int i = 0; i < N; ++i) {
		Temp.base[i] = log(x.base[i]);
	}
	return Temp;
}

R1 pow(double base, const R1 & exp) {
	R1 Temp(exp);
	int N = exp.getSize();
	for (int i = 0; i < N; ++i) {
		Temp.base[i] = pow(base, exp.base[i]);
	}
	return Temp;
}

R1 pow(const R1 & base, double exp) {
	R1 Temp(base);
	int N = base.getSize();
	for (int i = 0; i < N; ++i) {
		Temp.base[i] = pow(base.base[i], exp);
	}
	return Temp;
}

R1 pow(const R1 & base, const R1 & exp) {
	R1 Temp(base);
	if (base.getSize() == exp.getSize()) {
		int N = base.getSize();
		for (int i = 0; i < N; ++i) {
			Temp.base[i] = pow(base.base[i], exp.base[i]);
		}
		return Temp;
	}
	else {
		cout << "NOT THE SAME DIMENSIONS!" << endl;
		return Temp;
	}
}





R1 derivative(const R1& f, int Ver) {
	R1 Temp(f);
	int N = f.getSize();
	if (Ver == 1) { //Alg #1
		for (int i = 0; i < (N - 1); ++i) {
			Temp.base[i] = (f.base[i + 1] - f.base[i]) / (f.dx);
		}
		Temp.base[N - 1] = 0.75*Temp.base[N - 2];
		return Temp;
	}
	else if (Ver == 2) { //Alg #2
		for (int i = 1; i < (N - 1); ++i) {
			Temp.base[i] = (f.base[i + 1] - f.base[i - 1]) / (2 * f.dx);
		}
		Temp.base[0] = 0.25*Temp.base[1];
		Temp.base[N - 1] = 0.5*Temp.base[N - 2];
		return Temp;
	}
	else if (Ver == 3) { // Alg# 3
		for (int i = 2; i < (N - 2); ++i) {
			Temp.base[i] = 2 * (f.base[i + 1] - f.base[i - 1]) / (3 * f.dx) -
				(f.base[i + 2] - f.base[i - 2]) / (12 * f.dx);
		}
		Temp.base[1] = 0.5*Temp.base[2];
		Temp.base[0] = 0.25*Temp.base[1];
		Temp.base[N - 2] = 0.5*Temp.base[N - 3];
		Temp.base[N - 1] = 0.5*Temp.base[N - 2];
		return Temp;
	}
	else {
		cout << "INVALID ORDER!" << endl;
		return Temp;
	}
}

double integralD(double a, double b, const R1& f) {
	R1* Temp = new R1(f);
	int N = f.getSize();
	int N1 = Temp->getIndex(a);
	int N2 = Temp->getIndex(b);
	double Sum = 0.5*(f.base[N1] + f.base[N2])*f.dx;
	for (int i = N1 + 1; i < (N2 - 1); ++i) {
		Sum = Sum + f.base[i]*f.dx;
	}
	delete Temp;
	return Sum;
}

double int_dx(int i1, int i2, const R1 & f) {
	R1* Temp = new R1(f);
	int N = f.n;
	if (i2 > i1) {
		double Sum = 0.5*(f.base[i1] + f.base[i2])*f.dx;
		for (int i = i1 + 1; i < (i2 - 1); ++i) {
			Sum = Sum + f.base[i] * f.dx;
		}
		delete Temp;
		return Sum;
	}
	else {
		delete Temp;
		return 0.0;
	}
}

int R1::getIndex(double p) {
	if ((p > x1) && (p < x2)) {
		double N = (p - x1) / dx;
		int N2 = static_cast<int>(N);
		return N2;
	}
	else {
		cout << "NOT AN ELEMENT!" << endl;
		return 0;
	}
}


