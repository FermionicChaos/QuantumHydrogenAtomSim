//#include "stdafx.h"
#include "stdud.h"
#include "VecFieldR2.h"
#include "VecR2.h"
#include "R2.h"

VecFieldR2::VecFieldR2() {

}


VecFieldR2::~VecFieldR2() {
	//cout << "V field cleared" << endl;
}
//Wurks
VecFieldR2::VecFieldR2(const VecFieldR2 & vecf) {
	Vx = vecf.Vx; Vy = vecf.Vy;
}
//Wurks
VecFieldR2 & VecFieldR2::operator=(const VecFieldR2 & rhs) {
	this->Vx = rhs.Vx; this->Vy = rhs.Vy;
	return *this;
}
//Wurks
VecFieldR2::VecFieldR2(VecFieldR2 && vecf) {
	Vx = vecf.Vx; Vy = vecf.Vy;
}
//Wurks
VecFieldR2 & VecFieldR2::operator=(VecFieldR2 && rhs) {
	this->Vx = rhs.Vx; this->Vy = rhs.Vy;
	return *this;
}
//Wurks
VecFieldR2::VecFieldR2(R2& vx, R2& vy) {
	//WARNING! Make sure R2 arrays are same dimension and correctly bounded! No prevention in here yet!
	Vx = vx; Vy = vy;
}

VecR2 VecFieldR2::operator()(int i, int j) {
	return VecR2(Vx(i,j),Vy(i,j));
}

VecR2 VecFieldR2::operator()(double x, double y) {
	return VecR2(Vx(x,y),Vy(x,y));
}

VecR2 VecFieldR2::operator()(VecR2& pos) {
	return VecR2(Vx(pos),Vy(pos));
}

VecFieldR2 VecFieldR2::operator+(const VecFieldR2 & rhs) {
	return VecFieldR2(Vx + rhs.Vx, Vy + rhs.Vy);
}

VecFieldR2 VecFieldR2::operator-(const VecFieldR2 & rhs) {
	return VecFieldR2(Vx - rhs.Vx, Vy - rhs.Vy);
}

R2 VecFieldR2::operator*(const VecFieldR2 & rhs) {
	R2 temp;
	temp = Vx*rhs.Vx + Vy*rhs.Vy;
	return temp;
}

VecFieldR2 VecFieldR2::operator+(double rhs) {
	return VecFieldR2(Vx + rhs, Vy + rhs);
}
VecFieldR2 VecFieldR2::operator-(double rhs) {
	return VecFieldR2(Vx - rhs, Vy - rhs);
}
VecFieldR2 VecFieldR2::operator*(double rhs) {
	return VecFieldR2(Vx * rhs, Vy * rhs);
}
VecFieldR2 VecFieldR2::operator/(double rhs) {
	return VecFieldR2(Vx + rhs, Vy + rhs);
}
VecFieldR2 operator+(double lhs, const VecFieldR2 & vecf) {
	return VecFieldR2(lhs + vecf.Vx, lhs + vecf.Vy);
}
VecFieldR2 operator-(double lhs, const VecFieldR2 & vecf) {
	return VecFieldR2(lhs + vecf.Vx, lhs + vecf.Vy);
}
VecFieldR2 operator*(double lhs, const VecFieldR2 & vecf) {
	return VecFieldR2(lhs + vecf.Vx, lhs + vecf.Vy);
}
//Wurks
VecFieldR2 grad(const R2& f, int ver) {
	VecFieldR2 temp(d_dx(f, ver), d_dy(f, ver));
	return temp;
}
//Wurks
R2 div(const VecFieldR2 & vecf, int ver) {
	R2 temp(d_dx(vecf.Vx,ver) + d_dy(vecf.Vy,ver));
	return temp;
}
//Wurks
R2 curl(const VecFieldR2 & vecf, int ver) {
	R2 temp(d_dx(vecf.Vy, ver) - d_dy(vecf.Vx, ver));
	return temp;
}
