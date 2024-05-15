#include "Fluid.h"
#include <cmath>
#include <iostream>

Fluid::Fluid(int size_x, int size_y, float density, float h, float scale, float dt){
  this->size_x = size_x;
  this->size_y = size_y;
  this->count = this->size_x * this->size_y;
  this->density = density;
  this->h = h;
  this->dt = dt; 
  _h = 1.0 / h;
  _h2 = h / 2.0;
  sc = std::floor(scale * h);
  c.reserve(count);
  print();
  FluidCell test_c = FluidCell();
  setupCells();
}

void Fluid::extrapolate() {
  for (int i = 0; i < size_y; i++) {
    c[i * size_x + 0].v = c[i * size_x + 1].v;
    c[i * size_x + size_x - 1].v =
      c[i * size_x + size_x - 2].v;
  }
  for (int j = 0; j < size_x; j++) {
    c[0 * size_x + j].u = c[1 * size_x + j].u;
    c[(size_y - 1) * size_x + j].u =
      c[(size_y - 2) * size_x + j].u;
  }
}

void Fluid::projection() {
  float cp = (density * h) / dt;

  for (int i = 0; i < count; i++)
    c[i].p = 0;

  for (int iter = 0; iter < projection_iterations; iter++) {
    for (int i = 1; i < size_y - 1; i++) {
      for (int j = 1; j < size_x - 1; j++) {
        if (c[i * size_x + j].k == 0) continue;

        float divergence =
          c[(i + 1) * size_x + j].u -
          c[i * size_x + j].u +
          c[i * size_x + j + 1].v -
          c[i * size_x + j].v;

        int s_sum =
          c[(i - 1) * size_x + j].k +
          c[(i + 1) * size_x + j].k +
          c[i * size_x + j - 1].k +
          c[i * size_x + j + 1].k;

        float p = divergence /(float) s_sum;

        p *= overrelaxation;
        c[i * size_x + j].p += p * -cp;

        if (s_sum == 0) continue;
        c[i * size_x + j].u +=
          c[(i - 1) * size_x + j].k * p;

        c[i * size_x + j].v +=
          p * c[i * size_x + j - 1].k;

        c[(i + 1) * size_x + j].u -=
          p * c[(i + 1) * size_x + j].k;

        c[i * size_x + j + 1].v -=
          p * c[i * size_x + j + 1].k;
      }
    }
  }
}

void Fluid::advection() {
  for (int i = 0; i < count; i++) {
    c[i].u_new = c[i].u;
    c[i].v_new = c[i].v;
    c[i].s_new = c[i].s;
  }

  for (int i = 1; i < size_y - 1; i++) {
    for (int j = 1; j < size_x - 1; j++) {
      if (c[i * size_x + j].k == 0) continue;

      float x, y;
      // u component advection
      if (
        c[i * size_x + j].k == 1
      ) {
        float u = c[i * size_x + j].u;
        float v_avg =
          (c[i * size_x + j].v +
            c[i * size_x + j + 1].v +
            c[(i - 1) * size_x + j].v +
            c[(i - 1) * size_x + j + 1].v) /
          4;

        x = ((float) j) * h + _h2 - v_avg * dt;
        y = ((float) i) * h - u * dt;
        // calculate advection of velocities with backward euler
        c[i * size_x + j].u_new = sampleField(x, y, U);
      }

      // v component advection
      if (
        c[i * size_x + j].k == 1
      ) {
        float v = c[i * size_x + j].v;
        float u_avg =
          (c[i * size_x + j].u +
            c[i * size_x + j - 1].u +
            c[(i + 1) * size_x + j - 1].u +
            c[(i + 1) * size_x + j].u) /
          4;

        x = ((float) j) * h - v * dt;
        y = ((float) i) * h + _h2 - u_avg * dt;
        
        c[i * size_x + j].v_new = sampleField(x, y, V);
      }

      // smoke advection
      if (j < size_x - 1 && i < size_y - 1) {
        float u =
          (c[i * size_x + j].u +
            c[(i + 1) * size_x + j].u) *
          0.5;
        float v =
          (c[i * size_x + j].v +
            c[i * size_x + j + 1].v) *
          0.5;
        x = ((float) j) * h + _h2 - dt * v;
        y = ((float) i) * h + _h2 - dt * u; 
        c[i * size_x + j].s_new = sampleField(x, y, S);
      }
      // let vel = this.calcVelocity(i, j);
      // velSum += vel;
      // if (vel > this.currMaxVel) this.currMaxVel = vel;
      // if (!velSum) {
      //   console.log("error");
      //   let a;
      // }
    }
  }
  // this.avg_vel = velSum / ((size_x - 1) * (size_y - 1));
  for (int i = 0; i < count; i++) {
    c[i].u = c[i].u_new;
    c[i].v = c[i].v_new;
    c[i].s = c[i].s_new;
  }
}

float Fluid::sampleField(float x, float y, FieldType fieldType) {
  x = std::max(std::min(x, size_x * h), h);
  y = std::max(std::min(y, size_y * h), h);
  float dx = 0; // correction (0 or h/2)
  float dy = 0;
  if (fieldType == V) {
    dy = _h2;
  }
  if (fieldType == U) {
    dx = _h2;
  }
  if (fieldType == S) {
    dx = _h2;
    dy = _h2;
  }
  int x0 = std::min((int)std::floor((x - dx) * _h), size_x - 1);
  int y0 = std::min((int)std::floor((y - dy) * _h), size_y - 1);
  int x1 = std::min(x0 + 1, size_x - 1);
  int y1 = std::min(y0 + 1, size_y - 1);

  float w_x1 = (x - dx - (float)x0 * h) * _h;
  float w_x0 = 1 - w_x1;
  float w_y1 = (y - dy - (float)y0 * h) * _h;
  float w_y0 = 1 - w_y1;

  float sample_weighted = 0;
  if(fieldType == U) {
    sample_weighted = w_y0 * w_x0 * c[y0 * size_x + x0].u +
      w_y0 * w_x1 * c[y0 * size_x + x1].u +
      w_y1 * w_x0 * c[y1 * size_x + x0].u +
      w_y1 * w_x1 * c[y1 * size_x + x1].u;
  }
  if(fieldType == V) {
    sample_weighted = w_y0 * w_x0 * c[y0 * size_x + x0].v +
      w_y0 * w_x1 * c[y0 * size_x + x1].v +
      w_y1 * w_x0 * c[y1 * size_x + x0].v +
      w_y1 * w_x1 * c[y1 * size_x + x1].v;
  }
  if(fieldType == S) {
    sample_weighted = w_y0 * w_x0 * c[y0 * size_x + x0].s +
      w_y0 * w_x1 * c[y0 * size_x + x1].s +
      w_y1 * w_x0 * c[y1 * size_x + x0].s +
      w_y1 * w_x1 * c[y1 * size_x + x1].s;
  }
  // std::cout << sample_weighted << "\n";
  return sample_weighted;
}

void Fluid::get_pixel_array(int width, int height, sf::Uint8* pixels) {
  for (int i = 0; i < size_y; i++) {
    for (int j = 0; j < size_x; j++) {
      float p = c[i * size_x + j].p;
      float u = c[i * size_x + j].u;
      float v = c[i * size_x + j].v;
      float s = c[i * size_x + j].s;
      float k = c[i * size_x + j].k;
      int colors[3] = {255, 255, 255}; // r, g, b

      // if (scene.showPressure) colors = colorGradient(p, minP, maxP);
      if (k == 0) {
        colors[1] = 100;
        colors[2] = 200;
        colors[0] = 200;
      } else if (true) {
        colors[0] -= s * 255;
        colors[1] -= s * 255;
        colors[2] -= s * 255;

        if (colors[0] < 0) colors[0] = 0;
        if (colors[1] < 0) colors[1] = 0;
        if (colors[2] < 0) colors[2] = 0;
      }

      int x = j * sc; // in pixels
      int y = i * sc; // in pixels

      // process every pixel in a square at position (x,y), square size is h (meters) irl, h * scale on screen (pixels)
      for (int yp = y; yp < y + sc; yp++) {
        int index = (x + yp * width) * 4;
        for (int xp = x; xp < x + sc; xp++) {
          pixels[index++] = colors[0];
          pixels[index++] = colors[1];
          pixels[index++] = colors[2];
          pixels[index++] = 255;
        }
      }
    }
  }
}

void Fluid::setupCells() {
  std::cout << "Setting up cells\n";
  for (int i = 0; i < size_y; i++) {
    for (int j = 0; j < size_x; j++) {
      FluidCell new_cell = FluidCell();
      if (i == 0) {
        new_cell.k = 0;
      }
      if (i == size_y - 1) {
        new_cell.k = 0;
      }
      if (j == 0) {
        new_cell.k = 0;
      }
      if (j == size_x - 1) {
        new_cell.k = 1;
      }

      // obstacle
      if (j > size_x / 2 - 40 && j < size_x / 2 - 20 && i > size_y / 2 - 10 && i < size_y / 2 + 10) {
        new_cell.k = 0;
      }      
      // initial vel
      if (j == 1 || j == 0) {
        new_cell.v = initial_velocity;
			}
      // smoke
      if(j < 2 && i > size_y / 2 - 10 && i < size_y / 2 + 10) {
        new_cell.s = 1;
      }

      c.push_back(new_cell);
    }
  }
}

void Fluid::simulate() {
  // extrapolate();
  projection();
  advection();
}

void Fluid::print() {
  std::cout << "Printing data\n";
  std::cout << "Sim width: " << size_x * h << ", Sim height: " << size_y * h << "\n";
  std::cout << "Size x: " << size_x << ", Size y: " << size_y << "\n";
  std::cout << "Count: " << count << ", h: " << h << ", sc: " << sc << "\n";
  std::cout << "Density: " << density << ", dt: " << dt << "\n";
}

void Fluid::test() {
  c[4*size_x].print();
  c[4*size_x + 1].print();
  c[4*size_x + 10].print();
}

Fluid::~Fluid() {}