#include "point_mesh_distance.h"
#include "point_triangle_distance.h"
#include "igl/per_face_normals.h"


// idea
// iterate over the points in x
// iterate over the faces in FY
// run point to triangle distance, keeping only the smallest distance
// 	ensure that you keep track of (the distance, closest pt, and face of the
// 	smallest distance)
// add distance, position to closest pt to D and P respectively
// compute the normal and add it to M
void point_mesh_distance(
  const Eigen::MatrixXd & X,
  const Eigen::MatrixXd & VY,
  const Eigen::MatrixXi & FY,
  Eigen::VectorXd & D,
  Eigen::MatrixXd & P,
  Eigen::MatrixXd & N)
{
  Eigen::MatrixXd Normals;
  igl::per_face_normals(VY, FY, Normals);
  D.resize(X.rows());
  P.resize(X.rows(), 3);
  N.resize(X.rows(), 3);

  for (int i = 0; i < X.rows(); i++)
  {
    double shortest_d = 0;
    Eigen::RowVector3d position;
    int face_idx;

    for (int j = 0; j < FY.rows(); j++)
    {
      double d;
      Eigen::RowVector3d p;
      point_triangle_distance(X.row(i), VY.row(FY(j,0)), 
                              VY.row(FY(j, 1)), VY.row(FY(j, 2)), d, p);
      if (j==0)
      {
        shortest_d = d;
        position = p;
        face_idx = j;
      }
      else if (j != 0 && shortest_d > d)
      {
        shortest_d = d;
        position = p;
        face_idx = j;
      }
    }
    D(i) = shortest_d;
    P.row(i) = position;
    N.row(i) = Normals.row(face_idx);
  }
}
