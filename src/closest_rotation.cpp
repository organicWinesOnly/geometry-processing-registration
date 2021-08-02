#include "closest_rotation.h"
#include <Eigen/Dense>
#include <cassert>

void closest_rotation(
  const Eigen::Matrix3d & M,
  Eigen::Matrix3d & R)
{
  Eigen::JacobiSVD<Eigen::Matrix3d> M_svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
  Eigen::Matrix3d U = M_svd.matrixU();
  Eigen::Matrix3d V_T = M_svd.matrixV();  // returns the transpose of V
  int det_uv = (U * V_T).determinant();
  //assert(det_uv == 1 || det_uv == -1);

  Eigen::Matrix3d Q;
  Q = Eigen::Matrix3d::Identity();
  if (det_uv < 0)
  {
    Q(2, 2) = -1;
  }

  R = U * Q * V_T;
}
