#include <SFML/Graphics.hpp>
#include <cmath>
#include <iostream>
#include "Fluid.h"

int main()
{
  int width = 1200;
  int height = 800;
  float ratio = (float)height / (float)width;
  float h = 0.025;
  float sim_width = 4;
  float sim_height = sim_width * ratio;
  float scale = width / sim_width;
  float framerate = 60.0;
  float dt = 1 / framerate;

  Fluid fluid((int)std::floor(sim_width / h),(int) std::floor(sim_height / h), 1000, h, scale, dt);

  sf::RenderWindow window(sf::VideoMode(width, height), "Fluid Simulation");
  window.setFramerateLimit(framerate);

  while (window.isOpen())
  {
    sf::Event event;
    while (window.pollEvent(event))
    {
      if(event.type == sf::Event::Closed)
        window.close();
    }

    fluid.simulate();

    // drawing stuff
    window.clear(sf::Color::Black);
    sf::Image image;
    sf::Uint8* pixels = new sf::Uint8[width * height * 4];
    fluid.get_pixel_array(width, height, pixels);
    image.create(width, height, pixels);
    sf::Texture texture;
    texture.loadFromImage(image);
    
    sf::Sprite sprite(texture);
    window.draw(sprite);
    window.display();
  }

  return 0;
}