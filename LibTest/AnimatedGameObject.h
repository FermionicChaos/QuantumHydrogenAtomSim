#pragma once

#include "GameObject.h"
#include "R1.h"


using namespace std;
using namespace sf;

class AnimatedGameObject : public GameObject {
public:
	AnimatedGameObject();
	//~AnimatedGameObject();
	bool load(const string& filename);

	void draw(RenderWindow& window);

	void update(float deltaT);

	void setPosition(float x, float y);

	void setPosition(Vector2f R);

	void move(Vector2f v);

	Vector2f getPosition() const;

	void setIVec(float a, float b);

	float getHeight() const;

	float getWidth() const;

	void setScale(float scale);

private:
	Sprite m_sprite;
	Texture m_texture;
	string m_filename;
	bool m_valid = false;

	//float m; //Mass.
	//Vector2f p; //Momentum vector.
	Vector2f v; //Velocity vector.
	Vector2f R;
	Vector2f R0; //Initial position.
	float tt;
};

