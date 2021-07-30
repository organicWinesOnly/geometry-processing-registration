#include "point_to_plane_rigid_matching.h"
#include <iostream>
#include <cmath>
#include <Eigen/Geometry>
#include <Eigen/LU>

void axis_angle(
    const double theta,
    const Eigen::VectorXd w,
    Eigen::Matrix3d & R)
{
  Eigen::Matrix3d W;
  W << 0, -1 * w(2), w(1),
       w(2), 0, -1 * w(0),
       -1 * w(1), w(0), 0; 

  Eigen::Matrix3d W_sq = W * W;
  R = Eigen::Matrix3d::Identity() + sin(theta) * W + (1-cos(theta) ) * W_sq;
  
}

void point_to_plane_rigid_matching(
  const Eigen::MatrixXd & X,
  const Eigen::MatrixXd & P,
  const Eigen::MatrixXd & N,
  Eigen::Matrix3d & R,
  Eigen::RowVector3d & t)
{
  // create the A matrix and b vector
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(6, 6);
  Eigen::VectorXd b = Eigen::VectorXd::Zero(6);
  for (int i = 0; i < X.rows(); i++)
  {
    Eigen::Vector3d x = X.row(i);
    Eigen::Vector3d n = N.row(i);
    Eigen::VectorXd help_vector(6);
    help_vector.head(3)= x.cross(n);
    help_vector.tail(3)= N.row(i);
    b += help_vector * (N.row(i) * (P.row(i) - X.row(i) ).transpose() );
    for (int j = 0; j < 6; j++)
    {
      A.row(j) += help_vector(j) * help_vector.transpose();
    }
  }
  // solve Au = b for u in R^6
  Eigen::VectorXd u = A.inverse() * b;
  // the translation vector is given by t = (u[3], u[4], u[5])
  t = u.tail(3).transpose();
  // use the first 3 entries of u, the vector a, to compute the theta 
  Eigen::Vector3d a = u.head(3);
  // theta value and the w vector:
  //     theta = ||a||
  //     \hat{w} = a / ||a||
  double theta = a.norm();
  Eigen::Vector3d w = a / theta;
  // Then R will be given by axis_angle(theta, \hat{w})
  axis_angle(theta, w, R);
}
