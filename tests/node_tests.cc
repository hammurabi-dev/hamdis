// unit tests for Node class

#include <gtest/gtest.h>

#include <memory>
#include <hampixp.h>
#include <hampix.h>

const double cgs_pi = 3.14159265358979;

TEST(node, basic) {
  // default ctor
  Node<double> node_dft;
  EXPECT_EQ(node_dft.get_data(), double(0));
  EXPECT_EQ(node_dft.get_pointing().theta(), double(0));
  EXPECT_EQ(node_dft.get_pointing().phi(), double(0));
  EXPECT_EQ(node_dft.get_neighbor()->size(), int(0));
  
  // update data
  node_dft.update_data(0.2);
  EXPECT_EQ(node_dft.get_data(), double(0.2));
  
  // update pointing
  Hampixp new_point(1.2 * cgs_pi, 2.1 * cgs_pi);
  node_dft.update_pointing(new_point);
  EXPECT_NEAR(node_dft.get_pointing().theta(), double(0.2*cgs_pi), 1.0e-10);
  EXPECT_NEAR(node_dft.get_pointing().phi(), double(0.1*cgs_pi), 1.0e-10);
  
  // add neighbor
  auto node_neighbor = std::make_unique<Node<double>>();
  node_dft.add_neighbor(node_neighbor.get());
  EXPECT_EQ(node_dft.get_neighbor()->size(), int(1));
  EXPECT_EQ(node_dft.get_neighbor()->at(0)->get_data(), double(0));
  
  // move
  Node<double> node_mvconstr(std::move(node_dft));
  Node<double> node_mvassign;
  node_mvassign = std::move(node_mvconstr);
  
  // copy
  Node<double> node_cpconstr(node_mvassign);
  Node<double> node_cpassign;
  node_cpassign = node_cpconstr;
}
