#include "hausdorff_lower_bound.h"
#include "point_mesh_distance.h"
#include "random_points_on_mesh.h"

double hausdorff_lower_bound(
  const Eigen::MatrixXd & VX,
  const Eigen::MatrixXi & FX,
  const Eigen::MatrixXd & VY,
  const Eigen::MatrixXi & FY,
  const int n)
{
  // Build a set of points X
  Eigen::MatrixXd X(n, 3);

  random_points_on_mesh(n, VX, FY, X);
  // Find all the distances using point mesh distance function
  Eigen::VectorXd distance;
  Eigen::MatrixXd closest_pts;
  Eigen::MatrixXd norm_matrix;
  point_mesh_distance(X, VY, FY, distance, closest_pts, norm_matrix);

  return distance.maxCoeff();
}
