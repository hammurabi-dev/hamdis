// unit tests for Hampixp class

#include <cmath>
#include <gtest/gtest.h>
#include <hampixp.h>
#include <hamvec.h>
#include <iostream>
#include <memory>

const double cgs_pi = 3.14159265358979;

TEST(Hampixp, basic) {
  // init
  Hampixp ptr1;
  EXPECT_EQ(ptr1.Theta, double(0));
  EXPECT_EQ(ptr1.Phi, double(0));

  // return func
  EXPECT_EQ(ptr1.Theta, ptr1.theta());
  EXPECT_EQ(ptr1.Phi, ptr1.phi());

  // explicit init
  Hampixp ptr2(0.1, 0.2);
  EXPECT_EQ(ptr2.Theta, double(0.1));
  EXPECT_EQ(ptr2.Phi, double(0.2));

  // list init
  Hampixp ptr3 = {0.3, 0.4};
  EXPECT_EQ(ptr3.Theta, double(0.3));
  EXPECT_EQ(ptr3.Phi, double(0.4));

  // copy
  Hampixp ptr4(ptr2);
  EXPECT_EQ(ptr4.Theta, double(0.1));
  EXPECT_EQ(ptr4.Phi, double(0.2));
  ptr4 = ptr3;
  EXPECT_EQ(ptr4.Theta, double(0.3));
  EXPECT_EQ(ptr4.Phi, double(0.4));

  // norm
  Hampixp ptr5(1.2 * cgs_pi, 2.1 * cgs_pi);
  EXPECT_NEAR(ptr5.Theta, double(0.2 * cgs_pi), 1.0e-10);
  EXPECT_NEAR(ptr5.Phi, double(0.1 * cgs_pi), 1.0e-10);
  ptr5.theta(-1.2 * cgs_pi);
  EXPECT_NEAR(ptr5.Theta, double(0.8 * cgs_pi), 1.0e-10);
  ptr5.theta(-0.2 * cgs_pi);
  EXPECT_NEAR(ptr5.Theta, double(0.8 * cgs_pi), 1.0e-10);
  ptr5.theta(-0.9 * cgs_pi);
  EXPECT_NEAR(ptr5.Theta, double(0.1 * cgs_pi), 1.0e-10);
  ptr5.phi(-2.3 * cgs_pi);
  EXPECT_NEAR(ptr5.Phi, double(1.7 * cgs_pi), 1.0e-10);
  ptr5.phi(-0.8 * cgs_pi);
  EXPECT_NEAR(ptr5.Phi, double(1.2 * cgs_pi), 1.0e-10);
  ptr5.phi(-1.4 * cgs_pi);
  EXPECT_NEAR(ptr5.Phi, double(0.6 * cgs_pi), 1.0e-10);
}

TEST(Hampixp, vector) {
  Hampixp ptr1(0.3 * cgs_pi, 1.7 * cgs_pi);
  const auto v1 = ptr1.versor();
  EXPECT_NEAR(v1[0], std::sin(0.3 * cgs_pi) * std::cos(1.7 * cgs_pi), 1.0e-10);
  EXPECT_NEAR(v1[1], std::sin(0.3 * cgs_pi) * std::sin(1.7 * cgs_pi), 1.0e-10);
  EXPECT_NEAR(v1[2], std::cos(0.3 * cgs_pi), 1.0e-10);

  Hampixp ptr2(v1);
  EXPECT_NEAR(ptr2.Theta, double(0.3 * cgs_pi), 1.0e-10);
  EXPECT_NEAR(ptr2.Phi, double(1.7 * cgs_pi), 1.0e-10);
}
