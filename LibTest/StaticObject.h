#ifndef STATICOBJECT_H
#define STATICOBJECT_H
#pragma once

#include <SFML\Graphics.hpp>
#include "VecR2.h"
#include "GameObject.h"
using namespace std;
using namespace sf;


class StaticObject : public GameObject {
public:
	StaticObject();
	~StaticObject();

	bool load(const string& filename);

	void draw(RenderWindow& window);

	void update(float deltaT);

	void setPosition(float x, float y);

	void setPosition(Vector2f r);

	void move(Vector2f v);

	Vector2f getPosition() const;

	float getHeight() const;

	float getWidth() const;

	void setScale(float scale);

private:
	Sprite s_sprite;
	Texture s_texture;
	string s_filename;
	bool s_valid = false;
};
#endif // !STATICOBJECT_H


