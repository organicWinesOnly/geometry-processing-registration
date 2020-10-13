#include "point_mesh_distance.h"
#include "igl/per_face_normals.h"
#include <cmath>
#include <iostream>

void point_mesh_distance(
  const Eigen::MatrixXd & X,
  const Eigen::MatrixXd & VY,
  const Eigen::MatrixXi & FY,
  Eigen::VectorXd & D,
  Eigen::MatrixXd & P,
  Eigen::MatrixXd & N)
{
  N.resizeLike(X);
  P.resizeLike(X);
  D.resize(X.rows());
  Eigen::RowVector3d z;
  z.fill(1 / std::sqrt(3));
  Eigen::MatrixXd normals(FY.rows(), 3);

  // Find the normals of the every face in the mesh
  igl::per_face_normals(VY, FY, z, normals);

  // find the projected points
  Eigen::MatrixXd x_projections(normals.rows(), X.rows());
  x_projections = normals * X.transpose();

  // Find the closest face to each point in X
  for (int i = 0; i < X.rows(); i++)
  {
    int *idx;
    int index;

    idx = &index;
    double max_proj;
    max_proj = x_projections.col(i).maxCoeff(idx);
    N.row(i) = normals.row(*idx);
    P.row(i) =  max_proj * normals.row(*idx);
  }

  // distance calculation
  D = (X - P).rowwise().norm();
}
