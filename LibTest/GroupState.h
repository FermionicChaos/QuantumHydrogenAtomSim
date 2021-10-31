#ifndef GROUPSTATE_H
#define GROUPSTATE_H
#pragma once
#include "stdud.h"
#include <SFML/Graphics.hpp>
using namespace sf;

class GroupState {
private:
	unsigned char state;
	bool ds;

public:
	GroupState();
	~GroupState();
	//virtual GroupState(char a) = 0;

	virtual void load(int W, int H) = 0;
	virtual void Inp(Keyboard::Key, bool tf) = 0;
	virtual void update(float dt) = 0;
	virtual void draw(RenderWindow& window) = 0;

	virtual void clear() = 0;
	virtual char getState() = 0;
	virtual bool dS() = 0;
};
#endif // !GROUPSTATE_H

