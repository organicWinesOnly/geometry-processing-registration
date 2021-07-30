#include "point_triangle_distance.h"
#include "igl/per_face_normals.h"
#include <cassert>
#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <algorithm>


void point_triangle_distance(
  const Eigen::RowVector3d & x,
  const Eigen::RowVector3d & a,
  const Eigen::RowVector3d & b,
  const Eigen::RowVector3d & c,
  double & d,
  Eigen::RowVector3d & p)
{
  Eigen::RowVectorXd mean = (a + b+ c) / 3;
  Eigen::RowVector3d a_mean_center, b_mean_center, c_mean_center;
  Eigen::MatrixXd C(3,3);
  C.col(0) = (a-mean).transpose();
  C.col(1) = (b-mean).transpose();
  C.col(2) = (c-mean).transpose();
  Eigen::MatrixXd covarience;
  covarience = C.transpose() * C;

  Eigen::EigenSolver<Eigen::MatrixXd> eigensolve(covarience);
  Eigen::MatrixXd U(3, 2);
  U.col(0) = eigensolve.eigenvectors().col(0).real();
  U.col(1) = eigensolve.eigenvectors().col(1).real();

  // project edges and x onto the least square plane in 2D
   Eigen::MatrixXd planar_vertices;
   planar_vertices= U.transpose() * C;
   Eigen::VectorXd x_prime(2);
   x_prime = U.transpose() * x.transpose();
 
   // Use edge equation to find location of x_prime on least_sqaure plane
   // relative to the projected vertices
   double E12, E23, E31;
   E12 = (x_prime(0) -planar_vertices(0,0)) * (planar_vertices(1,1) - planar_vertices(1,0)) -
   	(x_prime(1) -planar_vertices(0,1)) * (planar_vertices(0,1) - planar_vertices(0,0));
 
   E23 = (x_prime(0) -planar_vertices(1,0)) * (planar_vertices(1,2) - planar_vertices(1,1)) -
     	 (x_prime(1) -planar_vertices(1,1)) * (planar_vertices(0,2) - planar_vertices(0,1));
 
   E31 = (x_prime(0) -planar_vertices(0,2)) * (planar_vertices(1,0) - planar_vertices(1,2)) -
   	 (x_prime(1) -planar_vertices(1,2)) * (planar_vertices(0,0) - planar_vertices(0,2));
 
 
  Eigen::VectorXd v12_double_prime = planar_vertices.col(1) - planar_vertices.col(0);
  Eigen::VectorXd q12_double_prime = x_prime - planar_vertices.col(0);
  Eigen::VectorXd v23_double_prime = planar_vertices.col(2) - planar_vertices.col(1);
  Eigen::VectorXd q23_double_prime = x_prime - planar_vertices.col(1);
  Eigen::VectorXd v31_double_prime = planar_vertices.col(0) - planar_vertices.col(2);
  Eigen::VectorXd q31_double_prime = x_prime - planar_vertices.col(2);

  double factor12 = std::sqrt(q12_double_prime.dot(v12_double_prime)) / (b-a).norm();
  Eigen::VectorXd pt_along_edge_12 =  factor12 * b + (1-factor12) * a;

  double factor23 = std::sqrt(q23_double_prime.dot(v23_double_prime)) / (c-b).norm();
  Eigen::VectorXd pt_along_edge_23 =  factor23 * c + (1-factor23) * b;

  double factor31 = std::sqrt(q31_double_prime.dot(v31_double_prime)) / (a-c).norm();
  Eigen::VectorXd pt_along_edge_31 =  factor31 * a + (1-factor31) * c;

  if (E12 < 0 and E23 < 0 and E31 < 0)
  {
    // x_prime lies in the projected triangle
    // take the mean of all the edge location
    p = 1/3 * (pt_along_edge_12 + pt_along_edge_23 + pt_along_edge_31).transpose();

  } else if (std::abs(E12) < 1e-15 or std::abs(E23) < 1e-15 or std::abs(E31) < 1e-15 )
  {
    // then x lies on one of the lines
    if (std::abs(E12) < 1e-15) 
    {
      p = pt_along_edge_12.transpose();
    } else if (std::abs(E23) < 1e-15) 
    {
      p = pt_along_edge_23.transpose();
    } else
    {
      p = pt_along_edge_31.transpose();
    }
  } else
  {
    // x lies outside of the triangle
    // take the minimum edge value, the pt lives closest to that edge
    if (E12 < std::min(E23, E31))
    {
      p = pt_along_edge_12.transpose();
    }else if (E23 < E31)
    {
      p = pt_along_edge_23.transpose();
    } else
    {
      p = pt_along_edge_31.transpose();
    }
  }
  d = (x-p).norm();
}

