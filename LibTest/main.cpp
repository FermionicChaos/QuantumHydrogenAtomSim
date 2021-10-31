#include "stdud.h"
#include "extra.h"

#include "Game.h"

using namespace std;

int main() {
	//Resolution.
	int hor = 512, ver = 512;
	int input;
	//cout << "Start Game?" << endl;
	//cin >> input;
	{
		Game game(hor, ver);
		//Game game();
		game.run();
	}
	//cout << "Close Window?" << endl;
	//cin >> input;
	return 0;
}

/*
cR2 x(1, X1, X2, RESX, Y1, Y2, RESY);
cR2 y(2, X1, X2, RESX, Y1, Y2, RESY);
cR2 z1, z2;
z1 = sin(pow(x, 4) + y*x)*tanh(x*y*y);
z2 = cos(pow(x, 4) + y*x)*jn(7, x*y*y);//*/

/*
CUDA seems to be about seven times faster than traditional computation. 
x card. 4096, y card. 4096, Computation on gpu finding these values was ~500 ms, on cpu, was ~3600 ms.
cudaR2v2 x(1, X1, X2, RESX, Y1, Y2, RESY);
cudaR2v2 y(2, X1, X2, RESX, Y1, Y2, RESY);
cout << "Domain Generated! " << clock.getElapsedTime().asMilliseconds() << "ms" << endl;
cin >> input;
clock.restart();
cudaR2v2 z = sqrt(x*x + y*y);
cout << "Computation Complete! " << clock.getElapsedTime().asMilliseconds() << "ms" << endl;
cin >> input;
clock.restart();
z.call2Host();
cout << "Data Called to host: " << clock.getElapsedTime().asMilliseconds() << " ms  " << z(4.0, 3.0) << endl;

cin >> input;
clock.restart();
R2 xx(1, X1, X2, RESX, Y1, Y2, RESY);
R2 yy(2, X1, X2, RESX, Y1, Y2, RESY);
cout << "Domain Generated! " << clock.getElapsedTime().asMilliseconds() << "ms" << endl;
cin >> input;
clock.restart();
R2 zz = sqrt(xx*xx + yy*yy);
cout << "Computation Complete! " << clock.getElapsedTime().asMilliseconds() << "ms" << endl;
cout << zz(4.0, 3.0) << endl;*/
