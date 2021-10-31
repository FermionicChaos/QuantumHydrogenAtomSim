#ifndef STAGE1_H
#define STAGE1_H
#pragma once



#include "GroupState.h"
#include <SFML/Graphics.hpp>
#include <SFML/OpenGL.hpp>
#include "GameObject.h"

#include "H_Atom_v2.h"
using namespace sf;

class Stage1 : public GroupState {
private:
	unsigned char state;
	bool ds;

	float buffer; //Time delay to distinguish between held and pulse.
	int hor, ver, N, pix;
	int select; //Select option.
	bool pause;

	H_Atom_v2 H2;

	GameObject *bg;
	//GameObject *fg;



	FloatRect *box; //Region of text objects.
	Text *txt; //Text objects.
	Font font; //Fonts used in game.
public:
	Stage1();
	~Stage1();
	Stage1(char a) { state = a; }

	void load(int W, int H);
	void Inp(Keyboard::Key key, bool tf);
	void update(float dt);
	void draw(RenderWindow& window);

	void clear();
	char getState();
	bool dS();
};
#endif // !STAGE1_H

