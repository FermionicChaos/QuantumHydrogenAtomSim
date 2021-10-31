#ifndef GeodesicRK4_H
#define GeodesicRK4_H

#pragma once

#include "R1.h"
#include "VecR4.h"

class GeodesicRK4 {
private:
	//Numerically solve the schwarzchild geodesic through RK4.
	R1 t, r, theta, phi; //Numerical functions storing the path of the object.
	R1 vt, vr, vtheta, vphi;
	VecR4 pos, vel; // Initial R4 position and velocity vectors.
	double dtau;
public:

	GeodesicRK4();
	~GeodesicRK4();
	//Constructor, allocate space for solving.
	GeodesicRK4(VecR4& Pos, VecR4& Vel, double dTau, int size);

	//Access Vector according to parameter.
	VecR4 GeodesicRK4::operator()(int I);
	VecR4 GeodesicRK4::operator()(double p0);

	//Default Schwarzschild motion solver RK4
	void solve();

	R1 getT() { return t; }
	R1 getX() { return r; }
	R1 getY() { return theta; }
	R1 getZ() { return phi; }
	VecR4 getR() { return pos; }
	VecR4 getV() { return vel; }

	void setT(R1& T) { t = T; }
	void setX(R1& X) { r = X; }
	void setY(R1& Y) { theta = Y; }
	void setZ(R1& Z) { phi = Z; }
	void setR(VecR4& R0) { pos = R0; }
	void setV(VecR4& V0) { vel = V0; }
};
#endif // !GeodesicRK4_H

