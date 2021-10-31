//
// Created by bswenson3 on 11/9/16.
//
#pragma once
#include "stdud.h"
#include "extra.h"

//#include <SFML/Graphics.hpp>
#include "GameObject.h"
#include "RegularGameObject.h"
#include "AnimatedGameObject.h"
#include "GroupState.h"

//#include "H_Atom_v2.h"

using namespace std;
using namespace sf;

class Game {
public:

    //Default is 640 by 480.
	Game();

	//Set resolution before start.
    Game(int W, int H);

    //Starts game.
    void run();

private:
	//Processes input.
    void processEvents();
    //update the game objects
    void update(Time deltaT);
    //draw the scene
    void render();
    //handle input from the user
    void handlePlayerInput(Keyboard::Key key, bool isDown);


    RenderWindow m_window;
	GroupState *g_state; //Game state group holder.
	//ContextSettings settings;

	unsigned char m_state; //Game State number (i.e, main menu, stage 1, etc...)
	bool dS, fullScr, running;
	int hor, ver;
	int input;
	/*
	0x00 Prologue
	0x01 Main Menu
	0x02 Stage 1
	*/

};


