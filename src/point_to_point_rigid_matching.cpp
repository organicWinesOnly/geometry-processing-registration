#include "point_to_point_rigid_matching.h"
#include "closest_rotation.h"

void point_to_point_rigid_matching(
  const Eigen::MatrixXd & X,
  const Eigen::MatrixXd & P,
  Eigen::Matrix3d & R,
  Eigen::RowVector3d & t)
{
  //R = Eigen::Matrix3d::Identity();
  //t = Eigen::RowVector3d::Zero();
  // find the centroids of P an X
  Eigen::Vector3d p_centroid;
  Eigen::Vector3d x_centroid;

  Eigen::VectorXd ones = Eigen::VectorXd::Ones(X.rows());
  double normalize_factor = ones.transpose() * ones;
  
  p_centroid = P.transpose() * ones;
  x_centroid = X.transpose() * ones;

  x_centroid = x_centroid / normalize_factor;
  p_centroid = p_centroid / normalize_factor;

  Eigen::MatrixXd P_mean = P;
  Eigen::MatrixXd X_mean = X;
  for (int i = 0; i < X.cols(); i++)
  {
    P_mean.col(i) = P_mean.col(i).array() - p_centroid(i); 
    X_mean.col(i) = X_mean.col(i).array() - x_centroid(i); 
  }

  Eigen::Matrix3d covarience_matrix = P_mean.transpose() * X_mean;
  closest_rotation(covarience_matrix, R);

  t = (p_centroid - R * x_centroid).transpose();
}

