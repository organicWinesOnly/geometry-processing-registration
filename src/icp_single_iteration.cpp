#include "icp_single_iteration.h"
#include "random_points_on_mesh.h"
#include "point_mesh_distance.h"
#include "point_to_point_rigid_matching.h"
#include "point_to_plane_rigid_matching.h"

void icp_single_iteration(
  const Eigen::MatrixXd & VX,
  const Eigen::MatrixXi & FX,
  const Eigen::MatrixXd & VY,
  const Eigen::MatrixXi & FY,
  const int num_samples,
  const ICPMethod method,
  Eigen::Matrix3d & R,
  Eigen::RowVector3d & t)
{
  // Build closest points matrix
  Eigen::MatrixXd P;
  // Build Normal matrix
  Eigen::MatrixXd N;
  // Build distance matrix
  Eigen::VectorXd D;
  // Build sample matrix
  Eigen::MatrixXd X;

  // Sample mesh
  random_points_on_mesh(num_samples, VX, FX, X);
  // Project onto mesh
  point_mesh_distance(X, VY, FY, D, P, N);
  // Update R and t
  if (method == ICP_METHOD_POINT_TO_POINT)
  {
    point_to_point_rigid_matching(X, P, R, t);
  } else
  {
    point_to_plane_rigid_matching(X, P, N, R, t);
  }
}
