#include "point_to_point_rigid_matching.h"
#include "closest_rotation.h"

void point_to_point_rigid_matching(
  const Eigen::MatrixXd & X,
  const Eigen::MatrixXd & P,
  Eigen::Matrix3d & R,
  Eigen::RowVector3d & t)
{
  // Compute centroid
  Eigen::RowVector3d p_centroid = P.colwise().sum() / P.rows();
  Eigen::RowVector3d x_centroid = X.colwise().sum() / X.rows();

  // Compute centoid distance matrix for X and P
  Eigen::MatrixXd p_bar(X.rows(), 3);
  Eigen::MatrixXd x_bar(X.rows(), 3);
  for (int i = 0; i < X.rows(); i++)
  {
    p_bar.row(i) = P.row(i) - p_centroid;
    x_bar.row(i) = X.row(i) - x_centroid;
  }
  
  // Copute covarience, M
  Eigen::Matrix3d M = p_bar.transpose() * x_bar;

  // Use closesnt relation function to compute the optimal R
  closest_rotation(M, R);

  // Calculate the optimal t
  t = p_centroid - (R * x_centroid.transpose()).transpose();
}

