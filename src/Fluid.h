#pragma once
#include <vector>
#include <cmath>
#include <SFML/Graphics.hpp>
#include <iostream>

class Fluid {
private:
  enum FieldType {U, V, S};

  struct FluidCell {
    float k = 0; // 0 - solid, 1 - fluid
    float u = 0; // vertical velocity
    float v = 0; // horizontal velocity
    float u_new = 0;
    float v_new = 0;
    float p = 0; // pressure
    float s = 0; // smoke amount
    float s_new = 0; 

    void print() {
      std::cout << "k: " << k << ", u: " << u << ", v: " << v << ", p: " << p << ", s: " << s << "\n";
    }
  };

  void extrapolate();
  void projection();
  void advection();
  float sampleField(float x, float y, FieldType field_type);
  void setupCells();
  
  float density;
  unsigned int size_x;
  unsigned int size_y;
  unsigned int count;
  float h;
  float _h;
  float _h2;
  float dt;
  int sc;
  int projection_iterations = 20;
  float initial_velocity = 2.0;
  float overrelaxation = 1.8;
  std::vector<FluidCell> c;

  
public:
  Fluid(unsigned int size_x, unsigned int size_y, float density, float h, float scale, float dt);
  ~Fluid();

  void simulate();
  void get_pixel_array(int width, int height, sf::Uint8* pixels);
  void print();
};