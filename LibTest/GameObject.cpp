//
// Created by bswenson3 on 11/9/16.
//

#include "GameObject.h"

using namespace std;
using namespace sf;

GameObject::GameObject() {

}
/*
bool GameObject::load(const string& filename) {
    if(m_texture.loadFromFile(filename)) {
        m_filename = filename;
        m_sprite.setTexture(m_texture);
        m_valid = true;
        return true;
    }
    return false;
}

void GameObject::draw(RenderWindow& window) {
    if(m_valid)
        window.draw(m_sprite);
}

void GameObject::move(Vector2f loc) {
    if(m_valid)
        m_sprite.move(loc);
}

void GameObject::setPosition(float x, float y) {
    if(m_valid)
        m_sprite.setPosition(x, y);
}

Vector2f GameObject::getPosition() const {
    if(m_valid)
        return m_sprite.getPosition();
    else
        return Vector2f(0, 0);
}

float GameObject::getHeight() const {
    return m_sprite.getLocalBounds().height;
}

float GameObject::getWidth() const {
    return m_sprite.getLocalBounds().width;
}

void GameObject::setScale(float scale) {
    if(m_valid)
        m_sprite.setScale(scale, scale);
}*/