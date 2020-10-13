#include "closest_rotation.h"
#include <Eigen/LU>
#include <Eigen/SVD>

void closest_rotation(
  const Eigen::Matrix3d & M,
  Eigen::Matrix3d & R)
{
  // Find the SVD of M
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
  Eigen::MatrixXd UV = svd.matrixU() * svd.matrixV().transpose();

  // compute the special matrix omega
  double omega_0;
  omega_0 = UV.determinant();
  Eigen::Matrix3d omega = Eigen::Matrix3d::Identity();
  omega(2, 2) = omega_0;

  // Solve for R
  R = svd.matrixU() * omega * svd.matrixV();
}
