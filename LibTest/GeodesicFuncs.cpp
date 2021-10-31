#include "GeodesicFuncs.h"

#include "stdud.h"
#include "R1.h"
#include "VecR4.h"

using namespace std;

double m = 1.0, M = 1.0, c = 1.0, G = 1.0;
//double m = 1.0, M = 5.972e24, c = 3e8, G = 6.67408e-11;
double E, L;

void EvalLE(VecR4& pos, VecR4& vel) {
	//Sets conserved quantities. E and L.
	L = m*pos.getX()*pos.getX()*vel.getZ();
	cout << "L = " << L << endl;
	double temp;
	//Uses normalization constraint to find dt/dTau. Feel free to jack off with velocity 3.
	temp = sqrt((1 - 2 * G*M / (pos.getX()*c*c))*(1+L*L/(m*m*c*c*pos.getX()*pos.getX())) + vel.getX()*vel.getX()/(c*c)) / (1-2*G*M/(pos.getX()*c*c));
	vel.setT(temp);
	cout << "dt/dtau = " << temp << endl;
	E = m*c*c*(1 - 2 * G*M / (pos.getX()*c*c))*vel.getT();
	cout << "E = " << E << endl;
}

double Rdot(int i, double r) {
	//Rearraged geodesic eq to numerically solve through rk4
	double rdot;
	double arg;
	arg = ((E*E - m*m*c*c*c*c) / (m*m*c*c)) + 2 * G*M / (r) - L*L / (m*m*r*r) + 2 * L*L*G*M / (m*m*c*c*r*r);
	rdot = sqrt(abs(arg));
	//cout << "Rdot = " << rdot << endl;
	return rdot;
}

double PHIdot(double r) {
	return (L/(m*r*r));
}

void PHI(const R1& r, const R1& phidot, R1& phi) {
	
	for (int i = 0; i < r.getSize(); ++i) {
		phi.getBase()[i] = (L/m)*int_dx(0, i, phidot) - phi.getBase()[0];
	}
}

double Tdot(double r) {
	double temp = (E/(m*c*c))*(r*c*c/(r*c*c - 2*G*M));
	return temp;
}

void T(R1 & tdot, R1 & t) {

	for (int i = 0; i < t.getSize(); ++i) {
		t.getBase()[i] = (L / m)*int_dx(0, i, tdot) - t.getBase()[0];
	}
}
//*/
