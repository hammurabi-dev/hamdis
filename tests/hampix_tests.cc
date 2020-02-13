// unit tests for Hampix class

#include <cmath>
#include <cstddef> // for std::size_t
#include <fstream>
#include <gtest/gtest.h>
#include <hampix.h>
#include <hampixp.h>
#include <iostream>
#include <memory>

const double cgs_pi = 3.14159265358979;

TEST(Hampix, basic) {
  // init
  Hampix<double> map_dft;
  EXPECT_EQ(map_dft.npix(), 0);

  // dft init with given size + copy + move + index + data
  Hampix<double> map_n(233);
  map_dft = map_n; // copy assign
  EXPECT_EQ(map_dft.npix(), 233);
  for (std::size_t i = 0; i < 233; ++i) {
    EXPECT_EQ(map_dft.index(i), i);
    EXPECT_EQ(map_dft.data(i), double(0.0));
    EXPECT_EQ(map_dft.pointing(i).phi(), double(0.0));
    EXPECT_EQ(map_dft.pointing(i).theta(), double(0.0));
  }
  Hampix<double> map_cpconstr(map_dft); // copy constr
  EXPECT_EQ(map_cpconstr.npix(), 233);
  for (std::size_t i = 0; i < 233; ++i) {
    EXPECT_EQ(map_cpconstr.index(i), i);
    EXPECT_EQ(map_cpconstr.data(i), double(0.0));
    EXPECT_EQ(map_cpconstr.pointing(i).phi(), double(0.0));
    EXPECT_EQ(map_cpconstr.pointing(i).theta(), double(0.0));
  }
  map_dft = std::move(map_cpconstr); // move assign
  EXPECT_EQ(map_dft.npix(), 233);
  for (std::size_t i = 0; i < 233; ++i) {
    EXPECT_EQ(map_dft.index(i), i);
    EXPECT_EQ(map_dft.data(i), double(0.0));
    EXPECT_EQ(map_dft.pointing(i).phi(), double(0.0));
    EXPECT_EQ(map_dft.pointing(i).theta(), double(0.0));
  }
  Hampix<double> map_mvconstr(std::move(map_n)); // move constr
  EXPECT_EQ(map_mvconstr.npix(), 233);
  for (std::size_t i = 0; i < 233; ++i) {
    EXPECT_EQ(map_mvconstr.index(i), i);
    EXPECT_EQ(map_mvconstr.data(i), double(0.0));
    EXPECT_EQ(map_mvconstr.pointing(i).phi(), double(0.0));
    EXPECT_EQ(map_mvconstr.pointing(i).theta(), double(0.0));
  }

  // pointing + index
  Hampixp point(0.3 * cgs_pi, 1.7 * cgs_pi);
  map_dft.pointing(23, point);
  map_dft.index(23, 1024);
  EXPECT_EQ(map_dft.index(23), std::size_t(1024));
  EXPECT_NEAR(map_dft.pointing(23).theta(), double(0.3 * cgs_pi), 1.0e-10);
  EXPECT_NEAR(map_dft.pointing(23).phi(), double(1.7 * cgs_pi), 1.0e-10);

  // constant init + print (requires visual inspection)
  Hampix<double> map_print(3, 0.2);
  for (std::size_t i = 0; i < 3; ++i) {
    EXPECT_EQ(map_print.data(i), double(0.2));
  }
  map_print.print();
}

TEST(Healmpix, basic) {
  // init

  Healmpix<double> map_dft;
  EXPECT_EQ(map_dft.npix(), 0);

  // constant init
  Healmpix<double> map_n(2);
  // nside + npix
  EXPECT_EQ(map_n.nside(), 2);
  EXPECT_EQ(map_n.npix(), 48);
  for (std::size_t i = 0; i < 48; ++i) {
    EXPECT_EQ(map_n.data(i), double(0));
    EXPECT_EQ(map_n.index(i), i);
  }

  // copy + move
  map_dft = map_n; // copy assign
  EXPECT_EQ(map_dft.nside(), 2);
  EXPECT_EQ(map_dft.npix(), 48);
  for (std::size_t i = 0; i < 48; ++i) {
    EXPECT_EQ(map_dft.data(i), double(0));
    EXPECT_EQ(map_dft.index(i), i);
  }
  Healmpix<double> map_cpconstr(map_dft); // copy constr
  EXPECT_EQ(map_cpconstr.nside(), 2);
  EXPECT_EQ(map_cpconstr.npix(), 48);
  for (std::size_t i = 0; i < 48; ++i) {
    EXPECT_EQ(map_cpconstr.data(i), double(0));
    EXPECT_EQ(map_cpconstr.index(i), i);
  }
  map_dft = std::move(map_cpconstr); // move assign
  EXPECT_EQ(map_dft.nside(), 2);
  EXPECT_EQ(map_dft.npix(), 48);
  for (std::size_t i = 0; i < 48; ++i) {
    EXPECT_EQ(map_dft.data(i), double(0));
    EXPECT_EQ(map_dft.index(i), i);
  }
  Healmpix<double> map_mvconstr(std::move(map_dft)); // move constr
  EXPECT_EQ(map_mvconstr.nside(), 2);
  EXPECT_EQ(map_mvconstr.npix(), 48);
  for (std::size_t i = 0; i < 48; ++i) {
    EXPECT_EQ(map_mvconstr.data(i), double(0));
    EXPECT_EQ(map_mvconstr.index(i), i);
  }

  // pointing calculation check with HEALPix Nside 2-32
  // reference pointing values are generated by healpy
  Healmpix<double> map_npower1(2);
  std::fstream infile1("reference/pointing_healpix_npower1.bin",
                       std::ios::in | std::ios::binary);
  if (infile1.is_open()) {
    for (std::size_t i = 0; i != map_npower1.npix(); ++i) {
      double tmp;
      // theta
      infile1.read(reinterpret_cast<char *>(&tmp), sizeof(double));
      EXPECT_NEAR(map_npower1.pointing(i).theta(), tmp, 1.0e-10);
      // phi
      infile1.read(reinterpret_cast<char *>(&tmp), sizeof(double));
      EXPECT_NEAR(map_npower1.pointing(i).phi(), tmp, 1.0e-10);
    }
    infile1.close();
  }
  Healmpix<double> map_npower2(4);
  std::fstream infile2("reference/pointing_healpix_npower2.bin",
                       std::ios::in | std::ios::binary);
  if (infile2.is_open()) {
    for (std::size_t i = 0; i != map_npower2.npix(); ++i) {
      double tmp;
      // theta
      infile2.read(reinterpret_cast<char *>(&tmp), sizeof(double));
      EXPECT_NEAR(map_npower2.pointing(i).theta(), tmp, 1.0e-10);
      // phi
      infile2.read(reinterpret_cast<char *>(&tmp), sizeof(double));
      EXPECT_NEAR(map_npower2.pointing(i).phi(), tmp, 1.0e-10);
    }
    infile2.close();
  }
  Healmpix<double> map_npower3(8);
  std::fstream infile3("reference/pointing_healpix_npower3.bin",
                       std::ios::in | std::ios::binary);
  if (infile3.is_open()) {
    for (std::size_t i = 0; i != map_npower3.npix(); ++i) {
      double tmp;
      // theta
      infile3.read(reinterpret_cast<char *>(&tmp), sizeof(double));
      EXPECT_NEAR(map_npower3.pointing(i).theta(), tmp, 1.0e-10);
      // phi
      infile3.read(reinterpret_cast<char *>(&tmp), sizeof(double));
      EXPECT_NEAR(map_npower3.pointing(i).phi(), tmp, 1.0e-10);
    }
    infile3.close();
  }
  Healmpix<double> map_npower4(16);
  std::fstream infile4("reference/pointing_healpix_npower4.bin",
                       std::ios::in | std::ios::binary);
  if (infile4.is_open()) {
    for (std::size_t i = 0; i != map_npower4.npix(); ++i) {
      double tmp;
      // theta
      infile4.read(reinterpret_cast<char *>(&tmp), sizeof(double));
      EXPECT_NEAR(map_npower4.pointing(i).theta(), tmp, 1.0e-10);
      // phi
      infile4.read(reinterpret_cast<char *>(&tmp), sizeof(double));
      EXPECT_NEAR(map_npower4.pointing(i).phi(), tmp, 1.0e-10);
    }
    infile4.close();
  }
  Healmpix<double> map_npower5(32);
  std::fstream infile5("reference/pointing_healpix_npower5.bin",
                       std::ios::in | std::ios::binary);
  if (infile5.is_open()) {
    for (std::size_t i = 0; i != map_npower5.npix(); ++i) {
      double tmp;
      // theta
      infile5.read(reinterpret_cast<char *>(&tmp), sizeof(double));
      EXPECT_NEAR(map_npower5.pointing(i).theta(), tmp, 1.0e-10);
      // phi
      infile5.read(reinterpret_cast<char *>(&tmp), sizeof(double));
      EXPECT_NEAR(map_npower5.pointing(i).phi(), tmp, 1.0e-10);
    }
    infile5.close();
  }
}