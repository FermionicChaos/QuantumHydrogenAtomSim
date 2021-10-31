#include "Menu.h"
//#include "GroupState.h"
#include <SFML/Graphics.hpp>
#include <iomanip>

#include "VecR2.h"

#include "GameObject.h"
#include "RegularGameObject.h"
#include "StaticObject.h"

using namespace std;
using namespace sf;

Menu::Menu() {
}

Menu::~Menu() {
	delete bg;

	delete[] txt;
	delete[] box;
}

void Menu::load(int W, int H) {
	hor = W; ver = H;
	cout << "Main Menu:" << endl;
	printf("Game State: %#x\n", state);

	N = 4;
	pix = ver/8;
	select = N;
	buffer = 0.0;
	ds = false;
	//state = 0x01;

	int x0, y0;
	x0 = hor / 2; y0 = ver / 2;

	if (!font.loadFromFile("images/fonts/AliciaWonderland.ttf")) {
		cout << "Error: Could not load fonts." << endl;
	}
	//Setup main text.
	title = Text();
	title.setFont(font);
	title.setString("Hydrogen Atom Simulation");
	title.setCharacterSize(ver/12);
	title.setColor(Color::Cyan);
	title.setStyle(Text::Italic);
	title.setOrigin(title.getLocalBounds().width/2, 0);
	title.setPosition(x0, 32);
	//Setup options.
	txt = new Text[N];
	box = new FloatRect[N];
	for (int i = 0; i < N; ++i) {
		txt[i].setFont(font);
	}
	txt[0].setString("Start");
	txt[1].setString("Options");
	txt[2].setString("Extra");
	txt[3].setString("Exit");
	for (int i = 0; i < N; ++i) {
		txt[i].setCharacterSize(pix);
		txt[i].setStyle(Text::Regular);
		box[i] = txt[i].getLocalBounds();
		txt[i].setOrigin(box[i].width/2, pix);
		txt[i].setPosition(x0,y0 + pix*(i+1) - (pix*N)/2);
	}

	bg = new StaticObject();
	bg->load("images/backgrounds/qm2.jpg");
	bg->setScale(0.6f);
}
//Region most operations will be done in each state.
void Menu::Inp(Keyboard::Key key, bool tf) {
	if ((tf)&&(buffer != 0.0)&&(buffer > 0.01)) {
		if (key == Keyboard::Up) {
			if ((select < N) && (select > 0)) {
				select -= 1;
			}
			else {
				//Bottom of the list.
				select = N - 1;
			}
		}
		if (key == Keyboard::Down) {
			if ((select < (N - 1)) && (select >= 0)) {
				select += 1;
			}
			else {
				//Top of the list.
				select = 0;
			}
		}
		if (key == Keyboard::Return) {
			if (select == 0) {
				state = 0x02;
				ds = true;
			}
			else if (select == 1) {

			}
			else if (select == 2) {

			}
			else if (select == 3) {
				state = 0xff;
				ds = true;
			}
			else {
				cout << "Inactive item" << endl;
			}
		}
	}
	else if ((tf) && (buffer != 0.0) && (buffer < 0.01)) {
		//Button held, buffer for better control.
	}
	else {
		buffer = 0.0;
	}
}

void Menu::update(float dt) {
	for (int i = 0; i < N; ++i) {
		if (i == select) {
			txt[i].setColor(Color::Blue);
		}
		else {
			txt[i].setColor(Color::Cyan);
		}
	}
	buffer += dt;
}

void Menu::draw(RenderWindow& window) {
	bg->draw(window);
	for (int i = 0; i < N; ++i) {
		window.draw(txt[i]);
	}
	window.draw(title);
}
//This will be used to specify methods that will delete pointer based objects
//Not cleared by default.
void Menu::clear() {

}

char Menu::getState() {
	return state;
}

bool Menu::dS() {
	return ds;
}
