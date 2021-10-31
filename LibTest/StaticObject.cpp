#include "StaticObject.h"



using namespace std;
using namespace sf;

StaticObject::StaticObject() {
}

StaticObject::~StaticObject() {
}


bool StaticObject::load(const string& filename) {
	if (s_texture.loadFromFile(filename)) {
		s_filename = filename;
		s_sprite.setTexture(s_texture);
		s_valid = true;
		return true;
	}
	return false;
}

void StaticObject::draw(RenderWindow& window) {
	if (s_valid)
		window.draw(s_sprite);
}

void StaticObject::update(float deltaT) {

}

void StaticObject::move(Vector2f loc) {
	if (s_valid)
		s_sprite.move(loc);
}

void StaticObject::setPosition(float x, float y) {
	if (s_valid)
		s_sprite.setPosition(x, y);
}

void StaticObject::setPosition(Vector2f r) {
	if (s_valid)
		s_sprite.setPosition(r.x, r.y);
}

Vector2f StaticObject::getPosition() const {
	if (s_valid)
		return s_sprite.getPosition();
	else
		return Vector2f(0, 0);
}

float StaticObject::getHeight() const {
	return s_sprite.getLocalBounds().height;
}

float StaticObject::getWidth() const {
	return s_sprite.getLocalBounds().width;
}

void StaticObject::setScale(float scale) {
	if (s_valid)
		s_sprite.setScale(scale, scale);
}
