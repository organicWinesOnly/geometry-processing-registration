#include "random_points_on_mesh.h"
#include "igl/cumsum.h"
#include "igl/doublearea.h"
#include <cassert>
#include <iostream>

void random_points_on_mesh(
  const int n,
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & X)
{
  X.resize(n, 3);

  // Build random vector
  Eigen::MatrixXd a_vec = Eigen::MatrixXd::Random(n, 2);
  // make all of a_vec values positive
  a_vec = a_vec.array().abs();
  //  make sure a1 + a2 <= 1
  Eigen::VectorXd a_vec_sum = a_vec.rowwise().sum();

  for (int i = 0; i < a_vec_sum.size(); i++)
  {
    if (a_vec_sum(i) > 1)
    {
      // replace a1 = 1 -a1, a2 = 1 - a2
      a_vec.row(i) = 1 - a_vec.row(i).array();
    }
  }

  // Area weighted random triangle sampling
  // Calculate the area of the triangle
  Eigen::VectorXd area(F.rows());
  igl::doublearea(V, F, area);
 
  double total_area = area.sum(); // total area of triangle mesh

  Eigen::VectorXd cum_sum(n);
  igl::cumsum(area / total_area, 1, cum_sum);
  Eigen::MatrixXd random_triangle(n, 9);
  
  int starting_pt = 0;
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < cum_sum.rows(); j++)
    {
      if (a_vec(i, 1) > cum_sum(j))
      {
	random_triangle.row(i) << V.row(F(j, 0)), V.row(F(j, 1)), V.row(F(j, 2));
	break;
      }
    }
  }

  // unifrom sampling of a single triangle
  Eigen::MatrixXd v1 = random_triangle.block(0, 0, n, 3);
  Eigen::MatrixXd v2 = random_triangle.block(0, 3, n, 3);
  Eigen::MatrixXd v3 = random_triangle.block(0, 6, n, 3);
  Eigen::MatrixXd a1(n, 3);
  Eigen::MatrixXd a2(n, 3);
  for (int i = 0; i < 3; i++)
  {
    a1.col(i) = a_vec.col(0);
    a2.col(i) = a_vec.col(1);
  }


  X = v1.array() + a1.array() * (v2 - v1).array() + a2.array() * (v3 - v1).array();
  // X = v1.array() + (v2 - v1).array();// + a_vec.col(1).array() * (v3 - v1).array();
}

