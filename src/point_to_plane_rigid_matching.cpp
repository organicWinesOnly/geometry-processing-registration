#include "point_to_plane_rigid_matching.h"
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Geometry>

void AxisAngle(
    Eigen::Vector3d omega,
    double theta,
    Eigen::Matrix3d & R)
{
  Eigen::Matrix3d W = Eigen::Matrix3d::Zero();
  W(0, 1) = - omega(2);
  W(0, 2) = omega(1);
  W(1, 0) = omega(2);
  W(1, 2) = - omega(0);
  W(2, 0) = - omega(1);
  W(2, 1) = omega(0);

  R = Eigen::Matrix3d::Identity() + std::sin(theta) * W + (1 - std::cos(theta)) * W;
}

void point_to_plane_rigid_matching(
  const Eigen::MatrixXd & X,
  const Eigen::MatrixXd & P,
  const Eigen::MatrixXd & N,
  Eigen::Matrix3d & R,
  Eigen::RowVector3d & t)
{
  // Build the LS problem
  // Build the A matrix
  Eigen::MatrixXd A(6, 6);
  Eigen::VectorXd b(6);
  A.setZero();
  b.setZero();

  for (int i = 0; i < X.rows(); i++)
  {
    Eigen::VectorXd temp_vec(6);
    Eigen::RowVector3d x_i = X.row(i);
    Eigen::RowVector3d n_i = N.row(i);

    temp_vec.head(3) = x_i.transpose().cross(n_i.transpose());
    temp_vec.tail(3) = n_i.transpose();
    Eigen::Vector3d diff = (P.row(i) - x_i).transpose();

    A += temp_vec * temp_vec.transpose();
    b += temp_vec * N.row(i) * diff;
  }

  // Solve least squares
  Eigen::VectorXd u(6);
  u = A.bdcSvd(Eigen::ComputeFullU| Eigen::ComputeFullV).solve(b);

  // Update parameters
  t = u.tail(3).transpose();
  double theta_new = u.head(3).norm();
  Eigen::Vector3d omega_new = u.head(3) / theta_new;
}
