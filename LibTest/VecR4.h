#ifndef VECR4_H
#define VECR4_H

#pragma once
#include <iostream>
class VecR4 {
private:
	double t, x, y, z;
public:
	VecR4();
	~VecR4();
	VecR4(VecR4& inp);
	VecR4& operator=(VecR4& rhs);
	VecR4(VecR4&& inp);
	VecR4& operator=(VecR4&& rhs);

	VecR4(double T, double X, double Y, double Z);

	friend std::ostream &operator<<(std::ostream &os, const VecR4& rhs);

	void setT(double T);
	void setX(double X);
	void setY(double Y);
	void setZ(double Z);

	double getT();
	double getX();
	double getY();
	double getZ();
};
#endif // !VECR4_H

