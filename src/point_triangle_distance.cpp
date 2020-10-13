#include "point_triangle_distance.h"
#include "igl/per_face_normals.h"
#include <iostream>

void point_triangle_distance(
  const Eigen::RowVector3d & x,
  const Eigen::RowVector3d & a,
  const Eigen::RowVector3d & b,
  const Eigen::RowVector3d & c,
  double & d,
  Eigen::RowVector3d & p)
{
  Eigen::Matrix3d vertex;
  vertex << a, b, c;
  Eigen::MatrixXi face(1, 3);
  face << 1, 2, 3;
  Eigen::RowVector3d normal;
  Eigen::RowVectorXd z;
  z << 1.0, 1.0, 1.0;
  igl::per_face_normals(vertex, face, z, normal);

  // Poject x onto the plane created by the triangle
  p = x.dot(normal) * normal; 

  // find distance
  Eigen::RowVector3d diff = x - p;
  d = diff.norm(); 
}
