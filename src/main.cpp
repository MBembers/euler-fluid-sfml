#include <SFML/Graphics.hpp>
#include <cmath>
#include <iostream>
#include "Fluid.h"

int main()
{
  int width = 1000;
  int height = 700;
  float ratio = (float)height / (float)width;
  float h = 0.02;
  float sim_width = 4;
  float sim_height = sim_width * ratio;
  float scale = width / sim_width;
  float framerate = 30.0;
  float dt = 1 / framerate;

  Fluid fluid(std::floor(sim_width / h), std::floor(sim_height / h), 1.0, h, scale, dt);

  // create the window
  sf::RenderWindow window(sf::VideoMode(width, height), "Fluid Simulation");
  window.setFramerateLimit(framerate);
  // create a clock to track the elapsed time
  sf::Clock clock;
  sf::Clock total_clock;
  total_clock.restart();

  // run the main loop
  while (window.isOpen())
  {
    // handle events
    sf::Event event;
    while (window.pollEvent(event))
    {
      if(event.type == sf::Event::Closed)
        window.close();
    }

    // make the particle system emitter follow the mouse
    // sf::Vector2i mouse = sf::Mouse::getPosition(window);

    // update clock
    sf::Time elapsed = clock.restart();

    fluid.simulate();
    // fluid.test();

    // drawing stuff
    window.clear(sf::Color::Green);
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