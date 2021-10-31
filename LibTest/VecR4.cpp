#include "VecR4.h"

#include <iostream>


VecR4::VecR4() {
	t = 0.0; x = 0.0;
	y = 0.0; z = 0.0;
}


VecR4::~VecR4() {

}

VecR4::VecR4(VecR4 & inp) {
	x = inp.getX(); y = inp.getY();
	t = inp.getT(); z = inp.getZ();
}

VecR4 & VecR4::operator=(VecR4 & rhs) {
	this->x = rhs.getX(); this->y = rhs.getY();
	this->t = rhs.getT(); this->z = rhs.getZ();
	return *this;
}

VecR4::VecR4(VecR4 && inp) {
	x = inp.getX(); y = inp.getY();
	t = inp.getT(); z = inp.getZ();
}

VecR4 & VecR4::operator=(VecR4 && rhs) {
	this->x = rhs.getX(); this->y = rhs.getY();
	this->t = rhs.getT(); this->z = rhs.getZ();
	return *this;
}

VecR4::VecR4(double T, double X, double Y, double Z) {
	t = T; x = X;
	y = Y; z = Z;
}

std::ostream &operator<<(std::ostream &os, VecR4 const &rhs) {
	return os << rhs.t << "*e_0 + " << rhs.x << "*e_1 + " << rhs.y << "*e_2 + " << rhs.z << "*e_3";
}

void VecR4::setT(double T) {
	t = T;
}

void VecR4::setX(double X) {
	x = X;
}

void VecR4::setY(double Y) {
	y = Y;
}

void VecR4::setZ(double Z) {
	z = Z;
}

double VecR4::getT() {
	return t;
}

double VecR4::getX() {
	return x;
}

double VecR4::getY() {
	return y;
}

double VecR4::getZ() {
	return z;
}
