// hampix is a multigrid spherical surface pixelization package
// based on HEALPix library
// in the very first version, we set only single layer multigrid
// while focusing on implementing basic functions
// but we leave the structure open to be upgraded
// use HEALPix nested ordering scheme only

#ifndef HAMMURABI_PIX_H
#define HAMMURABI_PIX_H

#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

// pixel info data type T
// minimal pixelization nside N
template <typename T> class Hampix {
public:
  // default constructor
  Hampix() = default;
  // copy constructor
  Hampix(const Hampix &m) { map.reset(new std::vector<Pix>(*(m.map.get()))); }
  // move constructor
  Hampix(Hampix &&m) noexcept : map{std::move(m.map)} {}
  // move assignment, no copy assignement since no const map
  Hampix &operator=(Hampix &m) noexcept {
    map = std::move(m.map);
    return *this;
  }
  // default destructor
  ~Hampix() = default;

#ifdef NDEBUG
private:
#endif
  // key structure of each pixel
  class Pix {
  public:
    // constructor
    Pix() : info{static_cast<T>(0)} {}
    ~Pix() = default;
    // pixel info
    T info;
    // pixelization HEALPix Nside tag
    unsigned int nside{0};
    // pixelization index at current Nside level
    std::size_t index{0};
    // pointer to nested parent pixel
    Pix *paxil{nullptr};
    // pointer to nested child pixels
    std::vector<Pix *> chxil{nullptr};
  };
  // map holder
  std::unique_ptr<std::vector<Pix>> map;
  // check HEALPix Nside setting
  bool nside_safety(const unsigned int &) const;
  // check HEALPix Npix setting
  bool npix_safety(const std::size_t &) const;
  // Nside to Npix calculator
  inline std::size_t nside2npix(const unsigned int &n) const {
    return static_cast<std::size_t>(12) * static_cast<std::size_t>(n) *
           static_cast<std::size_t>(n);
  }
  // Npix to Nside calculator
  inline unsigned int npix2nside(const std::size_t &n) const {
    return static_cast<unsigned int>(std::sqrt(n / 12));
  }

public:
  // initialize map with given Nside
  // all pix info assigned to zero
  void initialize(const unsigned int &);
  // initialize map from external std::vector
  void initialize(const std::vector<T> &);
  // print to content of each pix to screen
  void print() const;
  // pix2ang
};

template <typename T>
bool Hampix<T>::nside_safety(const unsigned int &suspect) const {
  const double spec{std::log2(suspect * double(suspect > 1))};
  return (spec == static_cast<unsigned int>(spec));
}

template <typename T>
bool Hampix<T>::npix_safety(const std::size_t &suspect) const {
  const double spec{std::log2(suspect * double(suspect > 1) / 12)};
  return (spec == static_cast<unsigned int>(spec));
}

template <typename T> void Hampix<T>::initialize(const unsigned int &N) {
  if (nside_safety(N)) { // Nside setting check
    map = std::make_unique<std::vector<Pix>>(12 * N * N);
    std::size_t idx{0};
    for (auto &i : *map) { // initialize nside and index
      i.nside = N;
      i.index = idx;
      idx++;
    }
  } else {
    throw std::runtime_error("wrong Nside setting");
  }
}

template <typename T> void Hampix<T>::initialize(const std::vector<T> &v) {
  const std::size_t length = v.size();
  if (npix_safety(length)) { // Npix setting check
    map = std::make_unique<std::vector<Pix>>(length);
    const unsigned int N = npix2nside(length);
    std::size_t idx{0};
    for (auto &i : *map) { // initialize nside and index
      i.nside = N;
      i.index = idx;
      i.info = v[idx];
      idx++;
    }
  } else {
    throw std::runtime_error("wrong Nside setting");
  }
}

template <typename T> void Hampix<T>::print() const {
  for (auto &i : *map) {
    std::cout << "info: " << i.info << "\t"
              << "Nside: " << i.nside << "\t"
              << "index: " << i.index << "\t" << std::endl;
  }
}

#endif
