#ifndef HAMMURABI_PIX_H
#define HAMMURABI_PIX_H

#include <cstddef> // for std::size_t
#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>
#include <hampixp.h>

// data type T
// key structure of each sample point
template <typename T>
class Node {
protected:
  // pointing
  Hampixp Pointing;
  // pixel data
  T Data {static_cast<T>(0)};
  // pixel index
  std::size_t Index=0;
  
public:
  // constructor
  Node() = default;
  // direct constr
  Node(const Hampixp &ptr, const T &val, const std::size_t &idx=0) {
    this->Pointing = ptr;
    this->Data = val;
    this->Index = idx;
  }
  // copy assign
  Node& operator=(const Node<T> &n) noexcept {
    this->Pointing = n.Pointing;
    this->Data = n.Data;
    this->Index = n.Index;
    return *this;
  }
  // move assign
  Node& operator=(Node<T> &&n) noexcept {
    this->Pointing = std::move(n.Pointing);
    this->Data = std::move(n.Data);
    this->Index = std::move(n.Index);
    return *this;
  }
  // copy constr
  Node(const Node<T> &n) noexcept {
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
  virtual Hampixp pointing() const {
    return this->Pointing;
  }
  // extract sky information
  virtual T data() const {
    return this->Data;
  }
  // extract sky index
  virtual std::size_t index() const {
    return this->Index;
  }
  // update sky position
  virtual void pointing(const Hampixp& new_pointing) {
    this->Pointing = new_pointing;
  }
  // update sky information
  virtual void data(const T& new_data) {
    this->Data = new_data;
  }
  // update ksy index
  virtual void index(const std::size_t& new_idx) {
    this->Index = new_idx;
  }
};

/*
// data type T
// hosts a list sampling points
template <typename T>
class Hampix {
protected:
  // map holder, unique pointer of vector of Node
  std::unique_ptr<std::vector<Node>> Map;
  
public:
  // dft constr
  Hampix() {
    thos->Map = std::make_unique<std::vector<Node>> ();
  }
  // copy constr
  Hampix(const Hampix &m) {
    this->Map.reset(new std::vector<Pix>(*(m.Map.get())));
  }
  // move constr
  Hampix(Hampix &&m) noexcept {
    this->Map = std::move(m.Map)
  }
  // move assign
  Hampix& operator=(Hampix &m) noexcept {
    this->Map = std::move(m.Map);
    return *this;
  }
  // copy assignment
  Hampix& operator=(const Hampix &&m)
  // dft destr
  ~Hampix() = default;
  
  // initialize map with given sample number
  // all Node info assigned to zero
  // no neighboring relation assigned
  void initialize(const unsigned int &);
  // initialize map with single value
  // all Node info assigned by the given value
  // no neighboring relation assigned
  void initialize(const T &);
  // print to content of each pix to screen
  void print() const;
  // pix2ang
};

template <typename T>
void Hampix<T>::initialize(const unsigned int &N) {
  map = std::make_unique<std::vector<Node>>(N);
}

template <typename T>
void Hampix<T>::initialize(const T &v) {
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

template <typename T>
void Hampix<T>::print() const {
  for (auto &i : *map) {
    std::cout << "info: " << i.info << "\t"
              << "Nside: " << i.nside << "\t"
              << "index: " << i.index << "\t" << std::endl;
  }
}
*/
#endif
