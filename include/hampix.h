// hampix is a pointing-based spherical surface pixelization package

#ifndef HAMMURABI_PIX_H
#define HAMMURABI_PIX_H

#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>
#include <hampixp.h>

// pixel data type T
// key structure of each sample point
template <typename T>
class Node {
protected:
  // pointing
  Hampixp point;
  // pixel data
  T data {static_cast<T>(0)};
  // unique pointer of vector of pointer to neighboring pixels
  std::unique_ptr<std::vector<Node<T> *>> neighbor;
  
public:
  // constructor
  Node() {
    neighbor = std::make_unique<std::vector<Node<T> *>>();
  }
  // direct constr
  Node(const Hampixp &ptr, const T &val) : point(ptr), data(val) {}
  // copy assign
  Node& operator=(const Node<T> &n) noexcept {
    point = n.point;
    data = n.data;
    neighbor.reset(new std::vector<Node<T> *>(n.ref_neighbor()));
    return *this;
  }
  // move assign
  Node& operator=(Node<T> &n) noexcept {
    point = std::move(n.point);
    data = std::move(n.data);
    neighbor = std::move(n.neighbor);
    return *this;
  }
  // copy constr
  Node(const Node<T> &n) noexcept {
    point = n.point;
    data = n.data;
    neighbor.reset(new std::vector<Node<T> *>(n.ref_neighbor()));
  }
  // move constr
  Node(Node<T> &&n) noexcept {
    point = std::move(n.point);
    data = std::move(n.data);
    neighbor = std::move(n.neighbor);
  }
  // destr
  virtual ~Node() = default;
  
  // extract sky position
  virtual Hampixp get_pointing() const{
    return point;
  }
  // extract sky information
  virtual T get_data() const{
    return data;
  }
  // extract neighbor ptr
  virtual std::vector<Node<T> *> * get_neighbor() const{
    return neighbor.get();
  }
  // extract neighbor content
  // designed for class copy semantics
  virtual std::vector<Node<T> *> ref_neighbor() const{
    return *neighbor;
  }
  // update sky position
  virtual void update_pointing(const Hampixp& new_point) {
    point = new_point;
  }
  // update sky information
  virtual void update_data(const T& new_data) {
    data = new_data;
  }
  // add neighbor
  virtual void add_neighbor(Node<T> *n) {
    neighbor->push_back(n);
  }
};

/*
// pixel info data type T
template <typename T>
class Hampix {
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
  // map holder, unique pointer of vector of Node
  std::unique_ptr<std::vector<Node>> map;
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
