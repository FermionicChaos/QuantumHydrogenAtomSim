#pragma once

#include "GameObject.h"

using namespace std;
using namespace sf;

class RegularGameObject : public GameObject {
public:
	RegularGameObject();
	//~RegularGameObject();

	bool load(const string& filename);

	void draw(RenderWindow& window);

	void update(float deltaT);

	void setPosition(float x, float y);

	void setPosition(Vector2f R);

	void move(Vector2f v);

	Vector2f getPosition() const;

	float getHeight() const;

	float getWidth() const;

	void setScale(float scale);

	void setR(bool inp) { m_right = inp; }
	void setL(bool inp) { m_left = inp; }
	void setU(bool inp) { m_up = inp; }
	void setD(bool inp) { m_down = inp; }
private:
	Sprite m_sprite;
	Texture m_texture;
	string m_filename;
	bool m_valid = false;

	float m_speed = 120.0f;
	bool m_left = false;
	bool m_right = false;
	bool m_up = false;
	bool m_down = false;

};

