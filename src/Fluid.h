#pragma once
#include <vector>
#include <cmath>
#include <SFML/Graphics.hpp>
#include <iostream>

class Fluid {
private:
  enum FieldType {U, V, S};

  struct FluidCell {
    float k; // 0 - solid, 1 - fluid
    float u; // vertical velocity
    float v; // horizontal velocity
    float u_new;
    float v_new;
    float p; // pressure
    float s; // smoke amount
    float s_new; 

    FluidCell() {
      k = 1;
      u = 0;
      v = 0;
      u_new = 0;
      v_new = 0;
      p = 0;
      s = 0;
      s_new = 0;
    }

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
  int size_x;
  int size_y;
  int count;
  float h;
  // inverse of h
  float _h;
  // half of h
  float _h2;
  float dt;
  int sc;
  int projection_iterations = 30;
  float initial_velocity = 2.0;
  float overrelaxation = 1.9;
  std::vector<FluidCell> c;

  
public:
  Fluid(unsigned int size_x, unsigned int size_y, float density, float h, float scale, float dt);
  ~Fluid();

  void simulate();
  void get_pixel_array(int width, int height, sf::Uint8* pixels);
  void print();
  void test();
};