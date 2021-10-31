#include "AnimatedGameObject.h"
#include "stdud.h"
#include "R1.h"

const double pi = 3.141592653589793238462643383279503;
const double e = 2.7182818284590452353602874713527;

AnimatedGameObject::AnimatedGameObject() {
	tt = 0;
}

//AnimatedGameObject::~AnimatedGameObject(){}

bool AnimatedGameObject::load(const string& filename) {
	if (m_texture.loadFromFile(filename)) {
		m_filename = filename;
		m_sprite.setTexture(m_texture);
		m_valid = true;
		return true;
	}
	return false;
}

void AnimatedGameObject::draw(RenderWindow& window) {
	if (m_valid)
		window.draw(m_sprite);
}

void AnimatedGameObject::update(float deltaT) {
	Vector2f R = getPosition();
	tt += deltaT;
	R.x = R0.x;
	R.y = 1000*gsl_sf_bessel_jl(5, tt) + R0.y;
	setPosition(R.x, R.y);
}


void AnimatedGameObject::setPosition(float x, float y) {
	if (m_valid)
		m_sprite.setPosition(x, y);
}

void AnimatedGameObject::setPosition(Vector2f R)
{
}

void AnimatedGameObject::move(Vector2f) {

}

Vector2f AnimatedGameObject::getPosition() const {
	if (m_valid)
		return m_sprite.getPosition();
	else
		return Vector2f(0, 0);
}

void AnimatedGameObject::setIVec(float a, float b) {
	R0.x = a; R0.y = b;
}

float AnimatedGameObject::getHeight() const {
	return m_sprite.getLocalBounds().height;
}

float AnimatedGameObject::getWidth() const {
	return m_sprite.getLocalBounds().width;
}

void AnimatedGameObject::setScale(float scale) {
	if (m_valid)
		m_sprite.setScale(scale, scale);
}

