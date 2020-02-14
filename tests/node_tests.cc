// unit tests for Node class

#include <cstddef> // for std::size_t
#include <gtest/gtest.h>
#include <hampix.h>
#include <hampixp.h>
#include <memory>

const double cgs_pi = 3.14159265358979;

TEST(node, basic) {
  // default ctor
  Node<double> node_dft;
  EXPECT_EQ(node_dft.data(), double(0));
  EXPECT_EQ(node_dft.pointing().theta(), double(0));
  EXPECT_EQ(node_dft.pointing().phi(), double(0));
  EXPECT_EQ(node_dft.index(), std::size_t(0));

  // update data
  node_dft.data(0.2);
  EXPECT_EQ(node_dft.data(), double(0.2));

  // update pointing
  Hampixp new_point(1.2 * cgs_pi, 2.1 * cgs_pi);
  node_dft.pointing(new_point);
  EXPECT_NEAR(node_dft.pointing().theta(), double(0.2 * cgs_pi), 1.0e-10);
  EXPECT_NEAR(node_dft.pointing().phi(), double(0.1 * cgs_pi), 1.0e-10);
  EXPECT_EQ(node_dft.index(), std::size_t(0));

  // update index
  node_dft.index(233);
  EXPECT_EQ(node_dft.index(), std::size_t(233));

  // move
  Node<double> node_mvconstr(std::move(node_dft));
  Node<double> node_mvassign;
  node_mvassign = std::move(node_mvconstr);
  EXPECT_NEAR(node_mvassign.pointing().theta(), double(0.2 * cgs_pi), 1.0e-10);
  EXPECT_NEAR(node_mvassign.pointing().phi(), double(0.1 * cgs_pi), 1.0e-10);
  EXPECT_EQ(node_mvassign.index(), std::size_t(233));

  // copy
  Node<double> node_cpconstr(node_mvassign);
  Node<double> node_cpassign;
  node_cpassign = node_cpconstr;
  EXPECT_NEAR(node_cpassign.pointing().theta(), double(0.2 * cgs_pi), 1.0e-10);
  EXPECT_NEAR(node_cpassign.pointing().phi(), double(0.1 * cgs_pi), 1.0e-10);
  EXPECT_EQ(node_cpassign.index(), std::size_t(233));
}
