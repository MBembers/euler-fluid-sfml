#include <SFML/Graphics.hpp>
#include "Fluid.h"

int main()
{
    float width = 800;
    float height = 500;
    // create the window
    sf::RenderWindow window(sf::VideoMode((int) width, (int) height), "Fluid Simulation");
    // create a clock to track the elapsed time
    sf::Clock clock;
    float ratio = height / width;
    float h = 0.1;
    float simWidth = 3;
    float simHeight = simWidth * ratio;
    float scale = width / simWidth;
    Fluid fluid(width / scale, height / scale, 1.0, 1.0, 0.1);

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

        // update it
        sf::Time elapsed = clock.restart();

        // draw it
        window.clear(sf::Color::White);
        window.display();
    }

    return 0;
}