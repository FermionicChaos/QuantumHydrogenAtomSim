#include "Stage1.h"
#include "StaticObject.h"
#include "cR2.h"

Stage1::Stage1() {

}


Stage1::~Stage1() {
	delete bg;
}

void Stage1::load(int W, int H) {
	hor = W; ver = H;
	cout << "Simulation:" << endl;
	printf("Game State: %#x\n", state);
	double X1 = -100.0, X2 = 100.0, Y1 = -100.0, Y2 = 100.0;
	int RESX = W, RESY = H;
	cR2 x(1, X1, X2, RESX, Y1, Y2, RESY);
	cR2 y(2, X1, X2, RESX, Y1, Y2, RESY);

	N = 7;
	ds = false;
	pause = false;
	//double* A = new double[N*(N + 1)*(2 * N + 1) / 6];
	//double* B = new double[N*(N + 1)*(2 * N + 1) / 6];

	H2.genIm(W, H);
	H2.genData(W, H, x, y, N);
	//H2.StateVec(A, B, x, y, N);

	//bg = new StaticObject();
	//bg->load("images/backgrounds/bg2.jpg");
	//bg->setScale(1.0f);
}

void Stage1::Inp(Keyboard::Key key, bool tf) {
}

void Stage1::update(float dt) {
		H2.compute(dt);
}

void Stage1::draw(RenderWindow & window) {
	H2.draw(window);
}

void Stage1::clear() {

}

char Stage1::getState() {
	return state;
}

bool Stage1::dS() {
	return ds;
}
