//
// Created by bswenson3 on 11/9/16.
//
#pragma once

#include <SFML/Graphics.hpp>
using namespace std;
using namespace sf;
/*
Generalized Game objects, divided into two classes: Static and Dynamic.

Static Objects: These objects remain stationary at all times, with the
caveat that they move with respect to the camera.

Dynamic Objects: These objects do not remain stationary and will react with the external
environment, along with other vector fields in the system.
*/
class GameObject {
public:
	//Defaut unused.
    GameObject();
	//load file from
    virtual bool load(const string& filename) = 0;

    virtual void draw(RenderWindow& window) = 0;

	virtual void update(float deltaT) = 0;

    virtual void setPosition(float x, float y) = 0;

	virtual void setPosition(Vector2f R) = 0;

    virtual void move(Vector2f v) = 0;

    virtual Vector2f getPosition() const = 0;

    virtual float getHeight() const = 0;

    virtual float getWidth() const = 0;

    virtual void setScale(float scale) = 0;

private:
    Sprite m_sprite;
    Texture m_texture;
    string m_filename;
    bool m_valid = false;
};


