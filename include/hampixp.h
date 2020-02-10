// reforged HEALPix pointing class
// spherical frame described by (theta,phi) in radian
// associated default Cartesian frame described by (x,y,z)
// and the default relation is
// theta = 0, z direction
// phi = 0, x direction
// phi = 90, y direction

#ifndef HAMMURABI_POINTING_H
#define HAMMURABI_POINTING_H

#include <cmath>
#include <hamvec.h>

const double cgspi = 3.14159265358979;

class Hampixp {
protected:
  double Theta{0};
  double Phi{0};
  // normalize the angles so that 0<=theta<=pi and 0<=phi<2*pi
  void norm();

public:
  Hampixp() = default;
  // copy constr
  Hampixp(const Hampixp &p) noexcept {
    this->Theta = p.theta();
    this->Phi = p.phi();
  }
  // copy assign
  Hampixp& operator=(const Hampixp &p) noexcept {
    this->Theta = p.theta();
    this->Phi = p.phi();
    return *this;
  }
  // move constr
  Hampixp(Hampixp &&p) noexcept {
    this->Theta = std::move(p.Theta);
    this->Phi = std::move(p.Phi);
  }
  // move assign
  Hampixp& operator=(Hampixp &&p) noexcept {
    this->Theta = std::move(p.Theta);
    this->Phi = std::move(p.Phi);
    return *this;
  }
  // explicit
  Hampixp(const double &t, const double &p) {
    this->Theta = t;
    this->Phi = p;
    norm();
  }
  // from Caretesian vector (not necessarily a versor)
  Hampixp(const Hamvec<3, double> &vec) {
    const double x{vec[0]};
    const double y{vec[1]};
    const double z{vec[2]};
    this->Theta = std::fmod(std::atan2(sqrt(x * x + y * y), z), cgspi);
    this->Theta += (this->Theta < 0) * cgspi;
    this->Phi = std::fmod(std::atan2(y, x), 2. * cgspi);
    this->Phi += (this->Phi < 0) * 2. * cgspi;
  }
  ~Hampixp() = default;

  // read in theta value
  void theta(const double &t) {
    this->Theta = std::fmod(t, cgspi);
    this->Theta += (this->Theta < 0) * cgspi;
  }
  // read in phi value
  void phi(const double &p) {
    this->Phi = std::fmod(p, 2. * cgspi);
    this->Phi += (this->Phi < 0) * 2. * cgspi;
  }
  // return theta value
  double theta() const { return this->Theta; }
  // return phi value
  double phi() const { return this->Phi; }
  // return to Cartesian versor
  Hamvec<3, double> versor() const;
};

void Hampixp::norm() {
  this->Theta = std::fmod(this->Theta, cgspi);
  this->Theta += (this->Theta < 0) * cgspi;
  this->Phi = std::fmod(this->Phi, 2. * cgspi);
  this->Phi += (this->Phi < 0) * 2. * cgspi;
}

Hamvec<3, double> Hampixp::versor() const {
  const double sth{std::sin(this->Theta)};
  return Hamvec<3, double>{sth * std::cos(this->Phi), sth * std::sin(this->Phi),
                           std::cos(this->Theta)};
}

#endif
