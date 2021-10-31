//
// Created by bswenson3 on 11/9/16.
//
#include "Game.h"

#include "stdud.h"
#include "extra.h"

#include "GroupState.h"
#include "Menu.h"
#include "Stage1.h"

using namespace std;
using namespace sf;

Game::Game() : m_window(VideoMode(640, 480), "Lame Game", Style::Close)
{
	//Game state xxxx
	m_state = 0x00; // Game State.
	dS = true; //State change (true/false)?
	fullScr = false; //Not set to full screen.
	hor = 640; // Horizontal resolution.
	ver = 480; //Vertical resolution.
	//Menu a;
	//g_state = new prologue(m_state);
}
/*
Try multithreading on calculations for speed up and removing the bottle, skip 
*/
Game::Game(int W, int H) : m_window(VideoMode(W, H), "H2 Quantum Simulator", Style::Close)
{
	/*
	ContextSettings settings;
	settings.depthBits = 24;
	settings.stencilBits = 8;
	settings.antialiasingLevel = 4;
	settings.majorVersion = 3;
	settings.minorVersion = 3;
	RenderWindow m_window(VideoMode(W, H), "Lame Game", Style::Close, settings);
	//*/


	//glEnable(GL_TEXTURE_2D);

	//Game state xxxx
	m_state = 0x01; // First Game State.
	dS = true; //State change (true/false)?
	fullScr = false; //Not set to full screen.
	hor = W; // Horizontal resolution.
	ver = H; //Vertical resolution.
	//g_state = new prologue(m_state);
}

void Game::run() {
	Clock clock;
	while (m_window.isOpen()) {

		//Change State.
		if (dS) {
			if (m_state == 0x00) {
				//cout << "Prologue:" << endl;
				//cout << hex << m_state << endl;

				//g_state = new prologue(m_state);
				//g_state->load(hor, ver);
			}
			else if (m_state == 0x01) {
				//delete g_state;
				g_state = new Menu(m_state);
				g_state->load(hor, ver);
			}
			else if (m_state == 0x02) {
				delete g_state;
				g_state = new Stage1(m_state);
				g_state->load(hor, ver);
			}
			else {
				//End game
				if (m_state != 0xff) {
					printf("State Transition Error Detected: %#x \n", m_state);
				}
				delete g_state;
				running = false;
				m_window.close();
			}
			//Reset transistion.
			dS = false;
		}
		else {
			Time deltaT = clock.restart();
			processEvents();
			if (dS)
				continue;
			update(deltaT);
			if (dS)
				continue;
			render();
		}
	}
}

void Game::processEvents() {
    Event event;
    while(m_window.pollEvent(event)) {
        switch(event.type) {
            case Event::KeyPressed:
                handlePlayerInput(event.key.code, true);
                break;
            case Event::KeyReleased:
                handlePlayerInput(event.key.code, false);
                break;
            case Event::Closed:
                m_window.close();
				running = false;
                break;
        }
    }
}

void Game::handlePlayerInput(Keyboard::Key key, bool isDown) {
	g_state->Inp(key, isDown);
}

//use time since last update to get smooth movement
void Game::update(Time deltaT) {
	g_state->update(deltaT.asSeconds());
}

void Game::render() {
    m_window.clear();

	g_state->draw(m_window);

	//glClear(GL_COLOR_BUFFER_BIT);

	m_window.display();
	//Check for change in state.
	m_state = g_state->getState();
	dS = g_state->dS();
}