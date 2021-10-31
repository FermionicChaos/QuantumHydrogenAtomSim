#include "RegularGameObject.h"

#include "GameObject.h"

using namespace std;
using namespace sf;

RegularGameObject::RegularGameObject() {

}


bool RegularGameObject::load(const string& filename) {
	if (m_texture.loadFromFile(filename)) {
		m_filename = filename;
		m_sprite.setTexture(m_texture);
		m_valid = true;
		return true;
	}
	return false;
}

void RegularGameObject::draw(RenderWindow& window) {
	if (m_valid)
		window.draw(m_sprite);
}

void RegularGameObject::update(float deltaT) {
	Vector2f movement(0.0f, 0.0f);
	if (m_up)
		movement.y -= m_speed;
	if (m_down)
		movement.y += m_speed;
	if (m_left)
		movement.x -= m_speed;
	if (m_right)
		movement.x += m_speed;

	move(movement * deltaT);
}

void RegularGameObject::move(Vector2f loc) {
	if (m_valid)
		m_sprite.move(loc);
}

void RegularGameObject::setPosition(float x, float y) {
	if (m_valid)
		m_sprite.setPosition(x, y);
}

void RegularGameObject::setPosition(Vector2f R) {
	if(m_valid)
		m_sprite.setPosition(R.x, R.y);
}

sf::Vector2f RegularGameObject::getPosition() const {
	if (m_valid)
		return m_sprite.getPosition();
	else
		return Vector2f(0, 0);
}

float RegularGameObject::getHeight() const {
	return m_sprite.getLocalBounds().height;
}

float RegularGameObject::getWidth() const {
	return m_sprite.getLocalBounds().width;
}

void RegularGameObject::setScale(float scale) {
	if (m_valid)
		m_sprite.setScale(scale, scale);
}