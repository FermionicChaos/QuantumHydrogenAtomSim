#ifndef MENU_H
#define MENU_H
#pragma once

#include "GroupState.h"
#include <SFML/Graphics.hpp>

#include "GameObject.h"
using namespace sf;

class Menu : public GroupState {
private:
	unsigned char state;
	bool ds; //State change.

	float buffer; //Time delay to distinguish between held and pulse.
	int hor, ver, N, pix;
	int select; //Select option.

	GameObject* bg;
	Text title; //Title
	FloatRect *box; //Region of text objects.
	Text *txt; //Text objects.
	Font font; //Fonts used in game.
public:
	Menu();
	~Menu();
	Menu(char a) { state = a; }

	void load(int W, int H);
	void Inp(Keyboard::Key key, bool tf);
	void update(float dt);
	void draw(RenderWindow& window);

	void clear();
	char getState();
	bool dS();
};
#endif // !MENU_H

