#ifndef VECFIELDR2_H
#define VECFIELDR2_H
#pragma once

#include "R2.h"
#include "VecR2.h"
class VecFieldR2 {
private:
	R2 Vx, Vy;
public:
	VecFieldR2();
	~VecFieldR2();
	VecFieldR2(const VecFieldR2& vecf);
	VecFieldR2& operator=(const VecFieldR2& rhs);
	VecFieldR2(VecFieldR2&& vecf);
	VecFieldR2& operator=(VecFieldR2&& rhs);

	VecFieldR2(R2& vx, R2& vy);

	VecR2 VecFieldR2::operator()(int i, int j);
	VecR2 VecFieldR2::operator()(double x, double y);
	VecR2 VecFieldR2::operator()(VecR2& pos);

	VecFieldR2 operator+(const VecFieldR2& rhs);
	VecFieldR2 operator-(const VecFieldR2& rhs);
	R2 VecFieldR2::operator*(const VecFieldR2& rhs);
	
	VecFieldR2 VecFieldR2::operator+(double rhs);
	VecFieldR2 VecFieldR2::operator-(double rhs);
	VecFieldR2 VecFieldR2::operator*(double rhs);
	VecFieldR2 VecFieldR2::operator/(double rhs);
	friend VecFieldR2 operator+(double lhs, const VecFieldR2& vecf);
	friend VecFieldR2 operator-(double lhs, const VecFieldR2& vecf);
	friend VecFieldR2 operator*(double lhs, const VecFieldR2& vecf);

	friend VecFieldR2 grad(const R2& f, int ver);
	friend R2 div(const VecFieldR2& vecf, int ver);
	friend R2 curl(const VecFieldR2& vecf, int ver);

	R2 getVX() { return Vx; }
	R2 getVY() { return Vy; }

	void setVX(R2& VX) { Vx = VX; }
	void setVY(R2& VY) { Vy = VY; }
};
#endif
