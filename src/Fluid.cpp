#include "Fluid.h"
#include <cmath>

Fluid::Fluid(unsigned int size_x, unsigned int size_y, float density, float h, float dt){
  this->size_x = size_x;
  this->size_y = size_y;
  this->count = size_x * size_y;
  this->density = density;
  this->h = h;
  this->_h = 1.0 / h;
  this->_h2 = h / 2.0;
  this->dt = dt; 
  c = std::vector<FluidCell>(count);
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

  for ( auto cell : c ) {
    cell.p = 0.0;
  }

  for (int iter = 0; iter < projection_iterations; iter++) {
    for (unsigned int i = 1; i < size_y - 1; i++) {
      for (unsigned int j = 1; j < size_x - 1; j++) {
        if (c[i * size_x + j].k == 0) continue;

        float divergence =
          c[(i + 1) *  + j].u -
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
  for ( auto cell : c ) {
    cell.u_new = cell.u;
    cell.v_new = cell.v;
    cell.s_new = cell.s;
  }

  for (size_t i = 1; i < size_y - 1; i++) {
    for (size_t j = 1; j < size_x - 1; j++) {
      if (c[i * size_x + j].k == 0) continue;

      float x, y;
      // u component advection
      if (
        c[(i - 1) * size_x + j].k != 0 &&
        j < size_x - 1
      ) {
        float u = c[i * size_x + j].u;
        float v_avg =
          (c[i * size_x + j].v +
            c[i * size_x + j + 1].v +
            c[(i - 1) * size_x + j].v +
            c[(i - 1) * size_x + j + 1].v) /
          4;

        x = j * h + _h2 - v_avg * dt;
        y = i * h - u * dt;
        // calculate advection of velocities with backward euler
        c[i * size_x + j].u_new = this->sampleField(
          x,
          y,
          U
        );
      }

      // v component advection
      if (
        c[i * size_x + j - 1].k != 0.0 &&
        i < size_y - 1
      ) {
        float v = c[i * size_x + j].v;
        float u_avg =
          (c[i * size_x + j].u +
            c[i * size_x + j - 1].u +
            c[(i + 1) * size_x + j - 1].u +
            c[(i + 1) * size_x + j].u) /
          4;

        x = j * h - v * dt;
        y = i * h + _h2 - u_avg * dt;
        c[i * size_x + j].v_new = this->sampleField(
          x,
          y,
          V
        );
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
        x = j * h + _h2 - dt * v;
        y = i * h + _h2 - dt * u;
        c[i * size_x + j].s_new = this->sampleField(
          x,
          y,
          S
        );
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
  for ( auto cell : c ) {
    cell.u = cell.u_new;
    cell.v = cell.v_new;
    cell.s = cell.s_new;
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
  float x0 = std::min(std::floor((x - dx) * _h), (float) size_x - 1);
  float y0 = std::min(std::floor((y - dy) * _h), (float) size_y - 1);
  float x1 = std::min(x0 + 1, (float) size_x - 1);
  float y1 = std::min(y0 + 1, (float) size_y - 1);

  float w_x1 = (x - dx - x0 * h) * _h;
  float w_x0 = 1 - w_x1;
  float w_y1 = (y - dy - y0 * h) * _h;
  float w_y0 = 1 - w_y1;

  if(fieldType == U) {
    return w_y0 * w_x0 * c[y0 * size_x + x0].u +
      w_y0 * w_x1 * c[y0 * size_x + x1].u +
      w_y1 * w_x0 * c[y1 * size_x + x0].u +
      w_y1 * w_x1 * c[y1 * size_x + x1].u;
  }
  if(fieldType == V) {
    return w_y0 * w_x0 * c[y0 * size_x + x0].v +
      w_y0 * w_x1 * c[y0 * size_x + x1].v +
      w_y1 * w_x0 * c[y1 * size_x + x0].v +
      w_y1 * w_x1 * c[y1 * size_x + x1].v;
  }
  else {
    return w_y0 * w_x0 * c[y0 * size_x + x0].s +
      w_y0 * w_x1 * c[y0 * size_x + x1].s +
      w_y1 * w_x0 * c[y1 * size_x + x0].s +
      w_y1 * w_x1 * c[y1 * size_x + x1].s;
  }
  // float sample_weighted =
  //   w_y0 * w_x0 * f[y0 * size_x + x0] +
  //   w_y0 * w_x1 * f[y0 * size_x + x1] +
  //   w_y1 * w_x0 * f[y1 * size_x + x0] +
  //   w_y1 * w_x1 * f[y1 * size_x + x1];
  // return sample_weighted;
}

void Fluid::setupCells() {
  for (size_t i = 0; i < count; i++) {
    for (size_t j = 0; j < count; j++) {
      c[i * size_x + j].k = 1.0;
      c[i * size_x + j].u = 0.0;
      c[i * size_x + j].v = 0.0;
      c[i * size_x + j].s = 0.0;

      if (i == 0) {
        c[i * size_x + j].k = 0;
      }
      if (i == size_y - 1) {
        c[i * size_x + j].k = 0;
      }
      if (j == 0) {
        c[i * size_x + j].k = 0;
      }
      if (j == size_x - 1) {
        c[i * size_x + j].k = 1;
      }

      if (j == 1) {
        c[i * size_x + j].v = initial_velocity;
        c[i * size_x + j].s = 1.0;
			}
    }
  }
}

void Fluid::simulate() {
  // this->extrapolate();
  this->projection();
  this->advection();
}

Fluid::~Fluid() {}