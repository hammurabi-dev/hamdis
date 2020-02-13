#ifndef HAMMURABI_PIX_H
#define HAMMURABI_PIX_H

#include <cassert>
#include <cmath>
#include <cstddef> // for std::size_t
#include <hampixp.h>
#include <iostream>
#include <memory>
#include <omp.h>
#include <stdexcept>
#include <vector>

const double halfpi = 1.57079632679490;

// data type T
// key structure of each sample point
template <typename T> class Node {
protected:
  // pointing
  Hampixp Pointing{0.0,0.0};
  // pixel data
  T Data{static_cast<T>(0)};
  // pixel index
  std::size_t Index{0};

public:
  // constructor
  Node() = default;
  // direct constr
  Node(const Hampixp &ptr, const T &val, const std::size_t &idx = 0) {
    this->Pointing = ptr;
    this->Data = val;
    this->Index = idx;
  }
  // copy assign
  Node &operator=(const Node<T> &n) {
    this->Pointing = n.Pointing;
    this->Data = n.Data;
    this->Index = n.Index;
    return *this;
  }
  // move assign
  Node &operator=(Node<T> &&n) noexcept {
    this->Pointing = std::move(n.Pointing);
    this->Data = std::move(n.Data);
    this->Index = std::move(n.Index);
    return *this;
  }
  // copy constr
  Node(const Node<T> &n) {
    this->Pointing = n.Pointing;
    this->Data = n.Data;
    this->Index = n.Index;
  }
  // move constr
  Node(Node<T> &&n) noexcept {
    this->Pointing = std::move(n.Pointing);
    this->Data = std::move(n.Data);
    this->Index = std::move(n.Index);
  }
  // destr
  virtual ~Node() = default;
  // extract sky position
  virtual Hampixp pointing() const { return this->Pointing; }
  // extract sky information
  virtual T data() const { return this->Data; }
  // extract sky index
  virtual std::size_t index() const { return this->Index; }
  // update sky position
  virtual void pointing(const Hampixp &new_pointing) {
    this->Pointing = new_pointing;
  }
  // update sky information
  virtual void data(const T &new_data) { this->Data = new_data; }
  // update ksy index
  virtual void index(const std::size_t &new_idx) { this->Index = new_idx; }
};

// data type T
// hosts a vector of Node objects
// the indices of Nodes are trivial
template <typename T> class Hampix {
protected:
  // map holder, unique pointer of vector of Nodes
  std::unique_ptr<std::vector<Node<T>>> Map;

public:
  // dft constr
  Hampix() { this->Map = std::make_unique<std::vector<Node<T>>>(); }
  // initialize map with given sample number N
  // Node's Data assigned by the given value
  // Node's Index assigned from 0 to N-1
  Hampix(const std::size_t &N, const T &v = static_cast<T>(0)) {
    this->Map = std::make_unique<std::vector<Node<T>>>(N);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (std::size_t i = 0; i < N; ++i) {
      this->Map->at(i).index(i);
      this->Map->at(i).data(v);
    }
  }
  // copy constr
  Hampix(const Hampix<T> &m) {
    this->Map.reset(new std::vector<Node<T>>(*(m.Map.get())));
  }
  // move constr
  Hampix(Hampix<T> &&m) noexcept { this->Map = std::move(m.Map); }
  // move assign
  Hampix &operator=(Hampix<T> &&m) noexcept {
    this->Map = std::move(m.Map);
    return *this;
  }
  // copy assignment
  Hampix &operator=(const Hampix<T> &m) {
    this->Map.reset(new std::vector<Node<T>>(*(m.Map.get())));
    return *this;
  }
  // dft destr
  virtual ~Hampix() = default;
  // extract node data
  virtual T data(const std::size_t &idx) const {
    return this->Map->at(idx).data();
  }
  // set node data
  virtual void data(const std::size_t &idx, const T &v) {
    this->Map->at(idx).data(v);
  }
  // extract node index
  virtual std::size_t index(const std::size_t &idx) const {
    return this->Map->at(idx).index();
  }
  // set node index
  virtual void index(const std::size_t &idx, const std::size_t &new_idx) {
    this->Map->at(idx).index(new_idx);
  }
  // extract node pointing
  virtual Hampixp pointing(const std::size_t &idx) const {
    return this->Map->at(idx).pointing();
  }
  // set node pointing
  virtual void pointing(const std::size_t &idx, const Hampixp &new_point) {
    this->Map->at(idx).pointing(new_point);
  }
  // extract map size
  virtual std::size_t npix() const { return this->Map->size(); }
  // print to content of each pix to screen
  virtual void print() const {
    std::cout << "... printing Hampix map information ..." << std::endl;
    // no need for multi-threading here
    for (auto i = this->Map->begin(); i < this->Map->end(); ++i) {
      std::cout << "index: " << i->index() << "\t"
                << "data: " << i->data() << "\t"
                << "pointing: theta " << i->pointing().theta() << " phi "
                << i->pointing().phi() << std::endl;
    }
  }
  // reset with given nside and clean up data
  virtual void reset(const std::size_t &n = 0) {
    // cleaning an used map with correct size
    if (n == 0 or this->Map->size() == n) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (std::size_t i = 0; i < this->Map->size(); ++i) {
        this->Map->at(i).data(0.0);
      }
    } else {
      // there 2 cases when Nside != n
      // 1, map to be recycled with wrong size
      // 2, empty map initialized by the default constr
      this->Map = std::make_unique<std::vector<Node<T>>>((n));
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (std::size_t i = 0; i < n; ++i) {
        this->Map->at(i).index(i);
        this->Map->at(i).pointing(Hampixp(0.0,0.0));
        this->Map->at(i).data(0.0);
      }
    }
  }
};

// using HEALPix RING ordering, with HEALPix Nside as power of 2
// also having trivial indexing
template <typename T> class Healmpix final : public Hampix<T> {
protected:
  // HEALPix nside_
  std::size_t Nside = 0;
  // HEALPix order_
  int Order = -1;
  // HEALPix npix_
  std::size_t Npix = 0;
  // HEALPix npface_
  std::size_t Npface = 0;
  // HEALPix ncap_
  std::size_t Ncap = 0;
  // HEALPix fact1_, fact2_
  double Fact1 = 0.0;
  double Fact2 = 0.0;
  // initialize constants
  // copied from healpix_base.cc SetNside function
  void prepare() {
    // order_
    this->Order = this->ilog2(this->Nside);
    // check order_
    if (this->Order < 0)
      throw std::runtime_error("unsupported HEALPix Nside");
    // npface_
    this->Npface = this->Nside * this->Nside;
    // npix_
    this->Npix = 12 * this->Npface;
    // ncap_
    this->Ncap = (this->Npface - this->Nside) << 1;
    // fact2
    this->Fact2 = 4. / this->Npix;
    // fact1
    this->Fact1 = (this->Nside << 1) * this->Fact2;
  }
  // copied from HEALPix cxxutils.h ilog2 function
  // Returns the largest integer n that fulfills 2^n<=arg.
  inline unsigned int ilog2(const std::size_t &n) {
    std::size_t tmp = n;
    unsigned int res = 0;
    while (tmp > 0x0000FFFF) {
      res += 16;
      tmp >>= 16;
    }
    if (tmp > 0x000000FF) {
      res |= 8;
      tmp >>= 8;
    }
    if (tmp > 0x0000000F) {
      res |= 4;
      tmp >>= 4;
    }
    if (tmp > 0x00000003) {
      res |= 2;
      tmp >>= 2;
    }
    if (tmp > 0x00000001) {
      res |= 1;
    }
    return res;
  }
  // copied from HEALPix cxxutils.h isqrt function
  // Returns the integer n, which fulfills n*n<=arg<(n+1)*(n+1).
  inline unsigned int isqrt(const std::size_t &n) {
    return unsigned(std::sqrt(static_cast<long double>(n) + 0.5));
  }
  // assigning correct pointing position
  // using HEALPix implementation, equivalent to HEALPix pix2ang function
  Hampixp fillpoint(const std::size_t &idx) {
    double z = 0.0;
    double phi = 0.0;
    double sth = 0.0;
    bool have_sth = false;
    // healpix RING ordering
    // copied from HEALPix healpix_base.cc pix2loc function
    // counted from North pole
    // North Polar cap
    if (idx < this->Ncap) {
      const std::size_t iring =
          (1 + static_cast<std::size_t>(this->isqrt(1 + 2 * idx))) >> 1;
      const std::size_t iphi = (idx + 1) - 2 * iring * (iring - 1);
      const double tmp = (iring * iring) * this->Fact2;
      z = 1.0 - tmp;
      if (z > 0.99) {
        sth = std::sqrt(tmp * (2.0 - tmp));
        have_sth = true;
      }
      phi = (iphi - 0.5) * halfpi / iring;
    }
    // Equatorial region
    else if (idx < (this->Npix - this->Ncap)) {
      const std::size_t nl4 = 4 * this->Nside;
      const std::size_t ip = idx - this->Ncap;
      const std::size_t tmp =
          (this->Order >= 0) ? ip >> (this->Order + 2) : ip / nl4;
      const std::size_t iring = tmp + this->Nside;
      const std::size_t iphi = ip - nl4 * tmp + 1;
      // 1 if iring+nside is odd, 1/2 otherwise
      const double fodd = ((iring + this->Nside) & 1) ? 1 : 0.5;
      // nside can be smaller than iring
      z = (static_cast<long double>(2 * this->Nside) -
           static_cast<long double>(iring)) *
          this->Fact1;
      phi = (iphi - fodd) * halfpi * 1.5 * this->Fact1;
    }
    // South Polar cap
    else {
      const std::size_t ip = this->Npix - idx;
      // counted from South pole
      const std::size_t iring =
          (1 + static_cast<std::size_t>(this->isqrt(2 * ip - 1))) >> 1;
      const std::size_t iphi = 4 * iring + 1 - (ip - 2 * iring * (iring - 1));
      const double tmp = (iring * iring) * this->Fact2;
      z = tmp - 1.0;
      if (z < -0.99) {
        sth = std::sqrt(tmp * (2.0 - tmp));
        have_sth = true;
      }
      phi = (iphi - 0.5) * halfpi / iring;
    }
    // copied from healpix_base.h pix2ang function
    return have_sth ? Hampixp(std::atan2(sth, z), phi)
                    : Hampixp(std::acos(z), phi);
  }

public:
  // dft constr
  Healmpix() : Hampix<T>() {}
  // initialize map with given HEALPix Nside
  // Node's Data assigned by the given value
  // Node's Index assigned from 0 to N-1
  Healmpix(const std::size_t &n, const T &v = static_cast<T>(0)) : Hampix<T>() {
    this->Nside = n;
    this->prepare();
    this->Map = std::make_unique<std::vector<Node<T>>>(
        static_cast<const std::size_t>(this->Npix));
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (std::size_t i = 0; i < this->Npix; ++i) {
      this->Map->at(i).index(i);
      this->Map->at(i).pointing(this->fillpoint(i));
      this->Map->at(i).data(v);
    }
  }
  // copy constr
  Healmpix(const Healmpix<T> &m) : Hampix<T>(m) {
    this->Nside = m.Nside;
    this->Order = m.Order;
    this->Npix = m.Npix;
    this->Npface = m.Npface;
    this->Ncap = m.Ncap;
    this->Fact1 = m.Fact1;
    this->Fact2 = m.Fact2;
  }
  // move constr
  Healmpix(Healmpix<T> &&m) : Hampix<T>(std::move(m)) {
    this->Nside = std::move(m.Nside);
    this->Order = std::move(m.Order);
    this->Npix = std::move(m.Npix);
    this->Npface = std::move(m.Npface);
    this->Ncap = std::move(m.Ncap);
    this->Fact1 = std::move(m.Fact1);
    this->Fact2 = std::move(m.Fact2);
  }
  // move assign
  Healmpix &operator=(Healmpix<T> &&m) noexcept {
    Hampix<T>::operator=(std::move(m));
    this->Nside = std::move(m.Nside);
    this->Order = std::move(m.Order);
    this->Npix = std::move(m.Npix);
    this->Npface = std::move(m.Npface);
    this->Ncap = std::move(m.Ncap);
    this->Fact1 = std::move(m.Fact1);
    this->Fact2 = std::move(m.Fact2);
    return *this;
  }
  // copy assignment
  Healmpix &operator=(const Healmpix<T> &m) {
    Hampix<T>::operator=(m);
    this->Nside = m.Nside;
    this->Order = m.Order;
    this->Npix = m.Npix;
    this->Npface = m.Npface;
    this->Ncap = m.Ncap;
    this->Fact1 = m.Fact1;
    this->Fact2 = m.Fact2;
    return *this;
  }
  // dft destr
  virtual ~Healmpix() = default;
  // nside
  std::size_t nside() const { return this->Nside; }
  // npix
  std::size_t npix() const override {
    assert(this->Npix == this->Map->size());
    return this->Npix;
  }
  // reset with given nside and clean up data
  void reset(const std::size_t &n = 0) override {
    // cleaning an used map with correct size
    if (n == 0 or this->Nside == n) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (std::size_t i = 0; i < this->Npix; ++i) {
        this->Map->at(i).data(0.0);
      }
    } else {
      // there 2 cases when Nside != n
      // 1, map to be recycled with wrong size
      // 2, empty map initialized by the default constr
      this->Nside = n;
      this->prepare();
      this->Map = std::make_unique<std::vector<Node<T>>>(
          static_cast<const std::size_t>(this->Npix));
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (std::size_t i = 0; i < this->Npix; ++i) {
        this->Map->at(i).index(i);
        this->Map->at(i).pointing(this->fillpoint(i));
        this->Map->at(i).data(0.0);
      }
    }
  }
};

#endif
