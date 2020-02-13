// hammurabi multi-dimensional (1D-3D) vector class based on std::vector

#ifndef HAMMURABI_VECTOR_H
#define HAMMURABI_VECTOR_H

#include <cassert>
#include <cmath>
#include <iostream>
#include <type_traits>
#include <vector>

template <int dim, typename T> class Hamvec {
protected:
  std::vector<T> Ele;

public:
  // default ctor
  // initialize with 0
  Hamvec() {
    switch (dim) {
    case 1:
      this->Ele = {T(0)};
      break;
    case 2:
      this->Ele = {T(0), T(0)};
      break;
    case 3:
      this->Ele = {T(0), T(0), T(0)};
      break;
    default:
      std::cerr << "unsupported dimension";
      break;
    }
  }
  virtual ~Hamvec() = default;
  // 1D vector
  Hamvec<dim, T>(const T &x) {
    assert(dim == 1);
    this->Ele.push_back(x);
  }
  // 2D vector
  Hamvec<dim, T>(const T &x, const T &y) {
    assert(dim == 2);
    this->Ele.push_back(x);
    this->Ele.push_back(y);
  }
  // 3D vector
  Hamvec<dim, T>(const T &x, const T &y, const T &z) {
    assert(dim == 3);
    this->Ele.push_back(x);
    this->Ele.push_back(y);
    this->Ele.push_back(z);
  }
  // copy ctor
  Hamvec<dim, T>(const Hamvec<dim, T> &v) { this->Ele = v.Ele; }
  // move ctor
  Hamvec<dim, T>(Hamvec<dim, T> &&v) { this->Ele = std::move(v.Ele); }
  // copy assign
  Hamvec<dim, T> &operator=(const Hamvec<dim, T> &v) noexcept {
    this->Ele = std::move(v.Ele);
    return *this;
  }
  // move assign
  Hamvec<dim, T> &operator=(Hamvec<dim, T> &&v) noexcept {
    this->Ele = std::move(v.Ele);
    return *this;
  }
  // constant operator []
  T operator[](const int &i) const { return this->Ele[i]; }
  // operator []
  T &operator[](const int &i) { return this->Ele[i]; }
  // get constant std::vector<T> Ele
  const std::vector<T> content() const { return this->Ele; }
  // get std::vector<T> Ele
  std::vector<T> &content() { return this->Ele; }
  // operator +
  // cast argument to the same template type
  template <typename R>
  Hamvec<dim, T> operator+(const Hamvec<dim, R> &v) const {
    Hamvec<dim, T> tmp(*this);
    for (unsigned int i = 0; i < dim; ++i) {
      tmp[i] += static_cast<T>(v[i]);
    }
    return tmp;
  }
  // operator +=
  // cast argument to the same template type
  template <typename R> Hamvec<dim, T> &operator+=(const Hamvec<dim, R> &v) {
    for (unsigned int i = 0; i < dim; ++i) {
      this->Ele[i] += static_cast<T>(v[i]);
    }
    return *this;
  }
  // operator -
  // cast argument to the same template type
  template <typename R>
  Hamvec<dim, T> operator-(const Hamvec<dim, R> &v) const {
    Hamvec<dim, T> tmp(*this);
    for (unsigned int i = 0; i < dim; ++i) {
      tmp[i] -= static_cast<T>(v[i]);
    }
    return tmp;
  }
  // operator -=
  // cast argument to the same template type
  template <typename R> Hamvec<dim, T> &operator-=(const Hamvec<dim, R> &v) {
    for (unsigned int i = 0; i < dim; ++i) {
      this->Ele[i] -= static_cast<T>(v[i]);
    }
    return *this;
  }
  // operator *
  // cast argument to the same template type
  template <typename R> Hamvec<dim, T> operator*(const R &s) const {
    Hamvec<dim, T> tmp(*this);
    for (unsigned int i = 0; i < dim; ++i) {
      tmp[i] *= static_cast<T>(s);
    }
    return tmp;
  }
  // operaotr *=
  // cast argument to the same template type
  template <typename R> Hamvec<dim, T> &operator*=(const R &s) {
    for (unsigned int i = 0; i < dim; ++i) {
      this->Ele[i] *= static_cast<T>(s);
    }
    return *this;
  }
  // operator /
  // cast argument to the same template type
  template <typename R> Hamvec<dim, T> operator/(const R &s) const {
    Hamvec<dim, T> tmp(*this);
    for (unsigned int i = 0; i < dim; ++i) {
      tmp[i] /= static_cast<T>(s);
    }
    return tmp;
  }
  // operator /=
  // cast argument to the same template type
  template <typename R> Hamvec<dim, T> &operator/=(const R &s) {
    assert(s != 0);
    for (unsigned int i = 0; i < dim; ++i) {
      this->Ele[i] /= static_cast<T>(s);
    }
    return *this;
  }
  // operator ==
  bool operator==(const Hamvec<dim, T> &v) {
    for (unsigned int i = 0; i < dim; ++i) {
      if (this->Ele[i] != v[i]) {
        return false;
      }
    }
    return true;
  }
  // operator !=
  bool operator!=(const Hamvec<dim, T> &v) {
    for (unsigned int i = 0; i < dim; ++i) {
      if (this->Ele[i] != v[i]) {
        return true;
      }
    }
    return false;
  }
  // vector length
  // cast into double type
  double length() const {
    double tmp{0};
    for (unsigned int i = 0; i < dim; ++i) {
      const double cache = static_cast<double>(this->Ele[i]);
      tmp += cache * cache;
    }
    return std::sqrt(tmp);
  }
  // vector squared length
  // cast into double type
  double lengthsq() const {
    double tmp{0};
    for (unsigned int i = 0; i < dim; ++i) {
      const double cache = static_cast<double>(this->Ele[i]);
      tmp += cache * cache;
    }
    return tmp;
  }
  // flip sign
  // should not be used to unsigned type
  void flip() {
    assert(std::is_signed<T>::value);
    for (unsigned int i = 0; i < dim; ++i) {
      this->Ele[i] *= static_cast<T>(-1.0);
    }
  }
  // versor
  // cast into double type
  Hamvec<dim, double> versor() const {
    Hamvec<dim, double> tmp;
    for (unsigned int i = 0; i < dim; ++i)
      tmp[i] = static_cast<double>(this->Ele[i]);
    const auto l2{tmp.lengthsq()};
    if (l2 != 0)
      tmp /= std::sqrt(l2);
    return tmp;
  }
  // inner product
  // cast into double type
  template <typename R> double dotprod(const Hamvec<dim, R> &v) const {
    double tmp{0};
    for (unsigned int i = 0; i < dim; ++i) {
      tmp += static_cast<double>(this->Ele[i]) * static_cast<double>(v[i]);
    }
    return tmp;
  }
  // cross product, works in 3D only
  // cast into double type
  template <typename R>
  Hamvec<dim, double> crossprod(const Hamvec<dim, R> &v) const {
    assert(dim == 3);
    Hamvec<dim, double> tmp;
    for (unsigned int i = 0; i < dim; ++i) {
      tmp[i] = (static_cast<double>(this->Ele[(i + 1) % 3]) *
                    static_cast<double>(v[(i + 2) % 3]) -
                static_cast<double>(this->Ele[(i + 2) % 3]) *
                    static_cast<double>(v[(i + 1) % 3]));
    }
    return tmp;
  }
  // osteam function
  friend std::ostream &operator<<(std::ostream &os, const Hamvec<dim, T> &v) {
    os << dim << "D vector: ";
    for (unsigned int i = 0; i < dim; ++i) {
      os << v[i] << "\t";
    }
    os << std::endl;
    return os;
  }
};

#endif
