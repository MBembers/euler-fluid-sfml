#pragma once
#include <vector>
#include <cmath>
#include <SFML/Graphics.hpp>

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
  int projection_iterations = 20;
  float initial_velocity = 2.0;
  float overrelaxation = 1.8;
  std::vector<FluidCell> c;

  
public:
  Fluid(unsigned int size_x, unsigned int size_y, float density, float h, float dt);
  ~Fluid();

  void simulate();
};