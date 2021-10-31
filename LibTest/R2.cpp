//#include "stdafx.h"
#include "stdud.h"
#include "VecR2.h"
#include "R2.h"
using namespace std;

const double pi = 3.141592653589793238;
const double e = 2.718281828459045235;
/*
int input;
cout << "ADDRESS: " << base << endl;
cin >> input;
*/
//Wurks
R2::R2() {
	cout << "DEFAULT" << endl;
	base = nullptr;
	double x1 = 0; double x2 = 0; double y1 = 0; double y2 = 0;
	double dx = 0; double dy = 0;
	int mem = 0; int m = 0; int n = 0;
}
//Wurks
R2::~R2() {
	cout << "ERASED" << endl;
	delete[] base;
}
//Wurks
// R2 ff(xx), or R2 ff = xx;
R2::R2(const R2& rhs) {
	cout << "COPY CONSTRUCT" << endl;
	dx = rhs.getDX(); dy = rhs.getDY();
	x1 = rhs.getX1(); x2 = rhs.getX2();
	y1 = rhs.getY1(); y2 = rhs.getY2();
	m = rhs.getIndex1(); n = rhs.getIndex2();
	base = new double[m*n];

	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			base[i + j*m] = rhs.getBase()[i + j*m];
		}
	}
}
//Wurks
// ff = xx
R2& R2::operator=(const R2& rhs) {
	cout << "COPY ASSIGN" << endl;
	if (base != rhs.getBase()) {
		this->dx = rhs.getDX(); this->dy = rhs.getDY();
		this->x1 = rhs.getX1(); this->x2 = rhs.getX2();
		this->y1 = rhs.getY1(); this->y2 = rhs.getY2();
		this->m = rhs.getIndex1(); this->n = rhs.getIndex2();
		delete[] this->base;
		base = new double[m*n];
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < n; ++j) {
				base[i + j*m] = rhs.getBase()[i + j*m];
			}
		}
		this->base = base;
		return *this;
	}
	else {
		return *this;
	}
}
//Wurks
// R2 ff(xx + yy), R2 ff = xx + yy
R2::R2(R2&& inp) {
	cout << "MOVE CONSTRUCT" << endl;
	base = inp.base;
	inp.base = nullptr;
	dx = inp.getDX(); dy = inp.getDY();
	x1 = inp.getX1(); x2 = inp.getX2();
	y1 = inp.getY1(); y2 = inp.getY2();
	m = inp.getIndex1(); n = inp.getIndex2();
}
// ff = xx + yy, ff = R2(xx + yy)
R2& R2::operator=(R2&& rhs) {
	delete[] base;
	base = rhs.getBase();
	dx = rhs.dx; dy = rhs.dy;
	x1 = rhs.x1; x2 = rhs.x2;
	y1 = rhs.y1; y2 = rhs.y2;
	m = rhs.m; n = rhs.n;
	rhs.base = nullptr;
	return *this;
}
//Wurks
R2::R2(int M, int N) {
	m = M; n = N;
	dx = 0.0; dy = 0.0;
	x1 = 0.0; x2 = 0.0; y1 = 0.0; y2 = 0.0;
	base = new double[M*N];
}
//Wurks
R2::R2(int k, double X1, double X2, int M, double Y1, double Y2, int N) {
	if ( (X2>X1) && (Y2>Y1) ) {
		m = M; n = N; x1 = X1; x2 = X2; y1 = Y1; y2 = Y2;
		base = new double[n*m];
		// FIX ROUND UP/DOWN CASTING
		double M2 = static_cast<double>(M);
		double N2 = static_cast<double>(N);
		dx = (x2 - x1) / M2;
		dy = (y2 - y1) / N2;
		if (k == 1) {
			double I = 0;
			for (int i = 0; i < m; ++i) {
				for (int j = 0; j < n; ++j) {
					base[i + j*m] = x1 + I*dx;
				}
				I += 1.0;
			}
			y2 = y1 + (n - 1)*dy; x2 = x1 + (m - 1)*dx;
		}
		else if (k == 2) {
			for (int i = 0; i < m; ++i) {
				double J = 0;
				for (int j = 0; j < n; ++j) {
					base[i + j*m] = y1 + J*dy;
					J += 1.0;
				}
			}
			y2 = y1 + (n - 1)*dy; x2 = x1 + (m - 1)*dx;
		}
		else {
			cout << "YOU GET GARBAGE, FUCK YOU!" << endl;
		}
	}
	else {
		m = M; n = N;
		dx = 0.0; dy = 0.0;
		x1 = 0.0; x2 = 0.0; y1 = 0.0; y2 = 0.0;
		base = new double[M*N];
	}
}

std::ostream &operator<<(std::ostream &os, R2 const &rhs) {
	return os << "x: " << "[" << rhs.getX1() <<
		", " << rhs.getX2() << "]" << " y: " << "[" <<
		rhs.getY1() << ", " << rhs.getY2() << "]" <<
		" Memory: " << rhs.getMEM() << " bytes";
}
//Wurks
double R2::operator()(int i, int j){
	if ((i < m) && (j < n) && (i >= 0) && (j >= 0)) {
		return base[i + j*m];
	}
	else {
		cout << "NOT AN ELEMENT IN THE SET!" << endl;
		return 0.0;
	}
}
//Wurks
double R2::operator()(double X, double Y) {
	if ((X > x1) && (X < x2) && (Y > y1) && (Y < y2)) {
		double Md = (X - x1) / dx; double Nd = (Y - y1) / dy;
		int M = static_cast<int>(Md); int N = static_cast<int>(Nd);
		return base[M + N*m];
	}
	else {
		cout << "NOT AN ELEMENT IN THE SET!" << endl;
		return 0.0;
	}
}

double R2::operator()(VecR2& pos) {
	if ((pos.getX() > x1) && (pos.getX() < x2) && (pos.getY() > y1) && (pos.getY() < y2)) {
		double Md = (pos.getX() - x1) / dx; double Nd = (pos.getY() - y1) / dy;
		int M = static_cast<int>(Md); int N = static_cast<int>(Nd);
		return base[M + N*m];
	}
	else {
		cout << "NOT AN ELEMENT IN THE SET!" << endl;
		return 0.0;
	}
}

//Next time, if conditions are met, use scoped object lol.
R2 R2::operator+(const R2& g) {
	cout << "Addition" << endl;
	if ( (m == g.getIndex1()) && (dx == g.getDX()) && (n == g.getIndex2()) && (dy == g.getDY()) ) {
		if ((x1 == g.getX1()) && (x2 == g.getX2()) && (y1 == g.getY1()) && (y2 == g.getY2())) {
			R2 temp(m, n);
			temp.setDX(dx); temp.setDY(dy);
			temp.setX1(x1); temp.setX2(x2);
			temp.setY1(y1); temp.setY2(y2);
			int N = m*n;
			double* Bout = temp.getBase();
			double* BL = base;
			double* BR = g.getBase();
			for (int i = 0; i < N; ++i) {
				Bout[i] = BL[i] + BR[i];
			}
			temp.setBase(Bout);
			return temp;
		}
		else {
			cout << "ERROR" << endl;
			R2 temp;
			return temp;
		}

	}
	else {
		cout << "ERROR" << endl;
		R2 temp;
		return temp;
	}
}

R2 R2::operator-(const R2& g) {
	cout << "Subtraction" << endl;
	if ((m == g.getIndex1()) && (dx == g.getDX()) && (n == g.getIndex2()) && (dy == g.getDY())) {
		if ((x1 == g.getX1()) && (x2 == g.getX2()) && (y1 == g.getY1()) && (y2 == g.getY2())) {
			R2 temp(m, n);
			temp.setDX(dx); temp.setDY(dy);
			temp.setX1(x1); temp.setX2(x2);
			temp.setY1(y1); temp.setY2(y2);
			int N = m*n;
			double* Bout = temp.getBase();
			double* BL = base;
			double* BR = g.getBase();
			for (int i = 0; i < N; ++i) {
				Bout[i] = BL[i] - BR[i];
			}
			temp.setBase(Bout);
			return temp;
		}
		else {
			cout << "ERROR" << endl;
			R2 temp;
			return temp;
		}

	}
	else {
		cout << "ERROR" << endl;
		R2 temp;
		return temp;
	}
}

R2 R2::operator*(const R2& g) {
	cout << "Multiplication" << endl;
	if ((m == g.getIndex1()) && (dx == g.getDX()) && (n == g.getIndex2()) && (dy == g.getDY())) {
		if ((x1 == g.getX1()) && (x2 == g.getX2()) && (y1 == g.getY1()) && (y2 == g.getY2())) {
			R2 temp(m, n);
			temp.setDX(dx); temp.setDY(dy);
			temp.setX1(x1); temp.setX2(x2);
			temp.setY1(y1); temp.setY2(y2);
			int N = m*n;
			double* Bout = temp.getBase();
			double* BL = base;
			double* BR = g.getBase();
			for (int i = 0; i < N; ++i) {
				Bout[i] = BL[i] * BR[i];
			}
			temp.setBase(Bout);
			return temp;
		}
		else {
			cout << "ERROR" << endl;
			R2 temp;
			return temp;
		}

	}
	else {
		cout << "ERROR" << endl;
		R2 temp;
		return temp;
	}
}

R2 R2::operator/(const R2& g) {
	cout << "Division" << endl;
	if ((m == g.getIndex1()) && (dx == g.getDX()) && (n == g.getIndex2()) && (dy == g.getDY())) {
		if ((x1 == g.getX1()) && (x2 == g.getX2()) && (y1 == g.getY1()) && (y2 == g.getY2())) {
			R2 temp(m, n);
			temp.setDX(dx); temp.setDY(dy);
			temp.setX1(x1); temp.setX2(x2);
			temp.setY1(y1); temp.setY2(y2);
			int N = m*n;
			double* Bout = temp.getBase();
			double* BL = base;
			double* BR = g.getBase();
			for (int i = 0; i < N; ++i) {
				if (BR[i] != 0.0) {
					Bout[i] = BL[i] / BR[i];
				}
				else {
					Bout[i] = 1.0;
					cout << "WARNING! Division by Zero attempted!" << endl;
				}
			}
			temp.setBase(Bout);
			return temp;
		}
		else {
			cout << "ERROR" << endl;
			R2 temp;
			return temp;
		}

	}
	else {
		cout << "ERROR" << endl;
		R2 temp;
		return temp;
	}
}

R2 R2::operator+(double rhs) {	
	R2 temp(*this);
	int N = temp.getSize();
	for (int i = 0; i < N; ++i) {
		temp.getBase()[i] = rhs + temp.getBase()[i];
	}
	return temp;
}

R2 R2::operator-(double rhs) {
	R2 temp(*this);
	int N = temp.getSize();
	for (int i = 0; i < N; ++i) {
		temp.getBase()[i] = temp.getBase()[i] - rhs;
	}
	return temp;
}

R2 R2::operator*(double rhs) {
	R2 temp(*this);
	int N = temp.getSize();
	for (int i = 0; i < N; ++i) {
		temp.getBase()[i] = temp.getBase()[i] * rhs;
	}
	return temp;
}

R2 R2::operator/(double rhs) {
	R2 temp(*this);
	if (rhs != 0) {
		int N = temp.getSize();
		for (int i = 0; i < N; ++i) {
			temp.getBase()[i] = temp.getBase()[i] / rhs;
		}
		return temp;
	}
	else {
		cout << "Division by zero not allowed!" << endl;
		return temp;
	}
}
// set [a b c d e ... n]
// a vs [b c d e ... n]
// then b vs [c d e ... n]
// (m*n-1)! operations
double R2::getMax() const {
//Return absolute max value.
	double temp, max = 0.0;
	for (int i = 0; i < m*n; ++i) {
		temp = base[i];
		if (abs(max) < abs(temp)) {
			max = base[i];
		}
	}
	return max;
}
//Wurks
bool R2::sameDim(const R2 & f, const R2 & g) {
	if ((f.m == g.m) && (f.dx == g.dx) && (f.n == g.n) && (f.dy == g.dy)) {
		if ((f.x1 == g.x1) && (f.x2 == g.x2) && (f.y1 == g.y1) && (f.y2 == g.y2)) {
			return true;
		}
		else {
			return false;
		}
	}
	else {
		return false;
	}
}


R2 operator+(double lhs, const R2& rhs) {
	R2 temp(rhs);
	int N = temp.getSize();
	for (int i = 0; i < N; ++i) {
		temp.getBase()[i] = lhs + temp.getBase()[i];
	}
	return temp;
}

R2 operator-(double lhs, const R2 & rhs) {
	R2 temp(rhs);
	int N = temp.getSize();
	for (int i = 0; i < N; ++i) {
		temp.getBase()[i] = lhs - temp.getBase()[i];
	}
	return temp;
}

R2 operator*(double lhs, const R2 & rhs) {
	R2 temp(rhs);
	int N = temp.getSize();
	for (int i = 0; i < N; ++i) {
		temp.getBase()[i] = lhs * temp.getBase()[i];
	}
	return temp;
}

R2 operator/(double lhs, const R2 & rhs) {
	R2 temp(rhs);
	int N = temp.getSize();
	for (int i = 0; i < N; ++i) {
		if (temp.getBase()[i] != 0) {
			temp.getBase()[i] = lhs / temp.getBase()[i];
		}
		else {
			temp.getBase()[i] = 1.0;
		}
	}
	return temp;
}

R2 sin(const R2& arg) {
	R2 temp(arg);
	int N = arg.getSize();
	for (int i = 0; i < N; ++i) {
		temp.getBase()[i] = sin(arg.getBase()[i]);
	}
	return temp;
}

R2 cos(const R2& arg) {
	R2 temp(arg);
	int N = arg.getSize();
	for (int i = 0; i < N; ++i) {
		temp.getBase()[i] = cos(arg.getBase()[i]);
	}
	return temp;
}

R2 tan(const R2& arg) {
	R2 temp(arg);
	int N = arg.getSize();
	for (int i = 0; i < N; ++i) {
		temp.getBase()[i] = tan(arg.getBase()[i]);
	}
	return temp;
}

R2 asin(const R2 & arg) {
	R2 temp(arg);
	int N = arg.getSize();
	for (int i = 0; i < N; ++i) {
		temp.getBase()[i] = asin(arg.getBase()[i]);
	}
	return temp;
}

R2 acos(const R2 & arg) {
	R2 temp(arg);
	int N = arg.getSize();
	for (int i = 0; i < N; ++i) {
		temp.getBase()[i] = acos(arg.getBase()[i]);
	}
	return temp;
}

R2 atan(const R2 & arg) {
	R2 temp(arg);
	int N = arg.getSize();
	for (int i = 0; i < N; ++i) {
		temp.getBase()[i] = atan(arg.getBase()[i]);
	}
	return temp;
}

R2 sinh(const R2& arg) {
	R2 temp(arg);
	int N = arg.getSize();
	for (int i = 0; i < N; ++i) {
		temp.getBase()[i] = sinh(arg.getBase()[i]);
	}
	return temp;
}

R2 cosh(const R2& arg) {
	R2 temp(arg);
	int N = arg.getSize();
	for (int i = 0; i < N; ++i) {
		temp.getBase()[i] = cosh(arg.getBase()[i]);
	}
	return temp;
}

R2 tanh(const R2& arg) {
	R2 temp(arg);
	int N = arg.getSize();
	for (int i = 0; i < N; ++i) {
		temp.getBase()[i] = tanh(arg.getBase()[i]);
	}
	return temp;
}

R2 d_dx(const R2 & f, int Ver) {
	R2 temp(f);
	int M = f.m, N = f.n;
	if (Ver == 1) {
		for (int j = 0; j < N; ++j) {
			for (int i = 0; i < M - 1; ++i) {
				temp.base[i + M*j] = (f.base[i+1 + M*j] - f.base[i + M*j])/(f.dx);
			}
			temp.base[M - 1 + M*j] = 0.75*temp.base[M - 1 + M*j];
		}
		return temp;
	}
	else if (Ver == 2) {
		for (int j = 0; j < N; ++j) {
			for (int i = 1; i < M - 1; ++i) {
				temp.base[i + M*j] = (f.base[i+1 + M*j] - f.base[i-1 + M*j]) / (2 * f.dx);
			}
			temp.base[M*j] = 0.75*temp.base[1 + M*j];
			temp.base[M - 1 + M*j] = 0.75*temp.base[M - 2 + M*j];
		}
		return temp;
	}
	else if (Ver == 3) {
		for (int j = 0; j < N; ++j) {
			for (int i = 2; i < M - 2; ++i) {
				temp.base[i + M*j] = 2 * (f.base[i+1 + M*j] - f.base[i-1 + M*j]) / (3 * f.dx)
					- (f.base[i+2 + M*j] - f.base[i-2 + M*j]) / (12 * f.dx);
			}
			temp.base[1 + M*j] = 0.75*temp.base[2 + M*j];
			temp.base[M*j] = 0.75*temp.base[1 + M*j];
			temp.base[M - 2 + M*j] = 0.75*temp.base[M - 3 + M*j];
			temp.base[M - 1 + M*j] = 0.75*temp.base[M - 2 + M*j];
		}
		return temp;
	}
	else {
		cout << "Invalid option" << endl;
		return temp;
	}
}

R2 d_dy(const R2 & f, int Ver) {
	R2 temp(f);
	int M = f.m, N = f.n;
	if (Ver == 1) {
		for (int i = 0; i < M; ++i) {
			for (int j = 0; j < N - 1; ++j) {
				temp.base[i + M*j] = (f.base[i + M*(j+1)] - f.base[i + M*j]) / (f.dy);
			}
			temp.base[i] = 0.75*temp.base[i + M];
		}
		return temp;
	}
	else if (Ver == 2) {
		for (int i = 0; i < M; ++i) {
			for (int j = 1; j < N - 1; ++j) {
				temp.base[i + M*j] = (f.base[i + M*(j+1)] - f.base[i + M*(j-1)]) / (2 * f.dx);
			}
			temp.base[i] = 0.75*temp.base[i + M];
			temp.base[i + M*(N-1)] = 0.75*temp.base[i + M*(N-2)];
		}
		return temp;
	}
	else if (Ver == 3) {
		for (int i = 0; i < N; ++i) {
			for (int j = 2; j < N - 2; ++j) {
				temp.base[i + M*j] = 2 * (f.base[i + M*(j+1)] - f.base[i + M*(j-1)]) / (3 * f.dx)
					- (f.base[i + M*(j+2)] - f.base[i + M*(j-2)]) / (12 * f.dx);
			}
			temp.base[i + M] = 0.75*temp.base[i + M*2];
			temp.base[i] = 0.75*temp.base[i + M];
			temp.base[i + M*(N-2)] = 0.75*temp.base[i + M*(N-3)];
			temp.base[i + M*(N-1)] = 0.75*temp.base[i + M*(N-2)];
		}
		return temp;
	}
	else {
		cout << "Invalid option" << endl;
		return temp;
	}
}

R2 normalize(R2& f) {
	R2 temp(f);
	double Max = f.getMax();
	temp = f / Max;
	return temp;
}
