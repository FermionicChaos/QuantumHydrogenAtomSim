#include "GeodesicRK4.h"

#include "stdud.h"

#include "GeodesuFuncs.h"

#include "R1.h"
#include "VecR4.h"

using namespace std;



GeodesicRK4::GeodesicRK4() {

}


GeodesicRK4::~GeodesicRK4() {

}

GeodesicRK4::GeodesicRK4(VecR4 & Pos, VecR4 & Vel, double dTau, int size) {
	pos = Pos; vel = Vel;
	t = R1(0.0, dTau*size, size); r = R1(0.0, dTau*size, size);
	theta = R1(0.0, dTau*size, size); phi = R1(0.0, dTau*size, size);
	vt = R1(0.0, dTau*size, size); vr = R1(0.0, dTau*size, size);
	vtheta = R1(0.0, dTau*size, size); vphi = R1(0.0, dTau*size, size);
	dtau = dTau;
}

VecR4 GeodesicRK4::operator()(int I) {
	return VecR4(t(I), r(I), theta(I), phi(I));
}

VecR4 GeodesicRK4::operator()(double p0) {
	return VecR4(t(p0), r(p0), theta(p0), phi(p0));
}

void GeodesicRK4::solve() {
	//Find R, and Rdot
	EvalLE(pos, vel);

	t.getBase()[0] = pos.getT();
	r.getBase()[0] = pos.getX();
	theta.getBase()[0] = pos.getY();
	phi.getBase()[0] = pos.getZ();

	vt.getBase()[0] = vel.getT();
	vr.getBase()[0] = vel.getX();
	vtheta.getBase()[0] = vel.getY();
	vphi.getBase()[0] = vel.getZ();

	int input;

	double k1, k2, k3, k4, dr;

	for (int i = 0; i < r.getSize() - 1; ++i) {
		if (i == 0) {
			k1 = vr.getBase()[0];
			//cout << k1 << " cocks " << Rdot(i, r(i)) << endl;
		}
		else {
			k1 = Rdot(i, r(i));
			//cout << "K1 = " << k1 << endl;
			if (vr.getBase()[i-1] > 0.0) {
				vr.getBase()[i] = k1;
			}
			else {
				vr.getBase()[i] = -k1;
			}
		}//*/
		//cout << k1 << endl;
		k2 = Rdot(i, r(i) + k1*dtau / 2);
		//cout << k2 << endl;
		k3 = Rdot(i, r(i) + k2*dtau / 2);
		k4 = Rdot(i, r(i) + k3*dtau);
		dr = dtau*(k1 + 2 * k2 + 2 * k3 + k4)/6;
		//cout << "abs dr = " << dr << endl;
		if (vr.getBase()[i] >= 0.0) {
			if ( (dr > 0.0001)&&(vr.getBase()[i] > 0.0001) ) {
				r.getBase()[i + 1] = r(i) + dr;
			}
			else {
				r.getBase()[i + 1] = r(i) - dr;
				vr.getBase()[i] = -k1;
			}
		}
		else if(vr.getBase()[i] < 0.0) {
				r.getBase()[i + 1] = r(i) - dr;
		}//*/
		else {
			continue;
		}
		//cout << " tau = " << i*dtau << endl;
		cout << "r = " << r.getBase()[i] << endl;
		cout << "dr/dtau = " << vr.getBase()[i] << endl;
		//cin >> input;
		if ((r(i) == 0.0) || (r(i) < 0.0)) {
			cout << "Terminated early, singularity reached." << endl;
			break;
		}
		//Do not forget the final element of rdot!
		//Complete R(tau).
	}
	//Generate rdot
	for (int i = 0; i < r.getSize(); ++i) {
		if (r(i) != 0.0) {
			vr.getBase()[i] = Rdot(i, r(i));
		}
		else {
			vr.getBase()[i] = 0.0;
		}
	}
	//Find phi, and phidot.
	for (int i = 0; i < r.getSize(); ++i) {
		vphi.getBase()[i] = PHIdot(r(i));
	}
	PHI(r, vphi, phi);

	for (int i = 0; i < r.getSize(); ++i) {
		vt.getBase()[i] = Tdot(r(i));
	}
	T(vt, t);

	theta = R1(0.0, 1.0, 10);
	theta = theta*0.0;
	vtheta = theta;


}
