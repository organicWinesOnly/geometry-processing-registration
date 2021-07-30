#include "random_points_on_mesh.h"
#include "igl/cumsum.h"
#include "igl/doublearea.h"
#include <cassert>
#include <iostream>
#include <random>

typedef std::default_random_engine G;
typedef std::uniform_real_distribution<double> Real_D;
// Idea
// create 2 functions
//   (1) for randomly sampling a triangle
//   (2) for randomly sampling the triangle indicies


// pt_on_triangle Returns a point in a triangle
//  input: vertices -> 3x3 matrix containing vertices of a triangle
//  output: position (x,y,z)
//
//  randomly generate 2 numbers from [0,1] and ensure the sum is less than 1
//  use the formula to find the position
void pt_on_triangle(
    const Eigen::Vector3d v1, 
    const Eigen::Vector3d v2, 
    const Eigen::Vector3d v3, 
    Eigen::Vector3d & position);


// triangle_in_mesh Compute binary search over triangle mesh areas
//  input: cumsum -> cumulative sum of triangle mesh areas
//         a -> random number from (0,1) 
//         start, stop -> beginning and end of array
//
//  output: triangle mesh index
int triangle_in_mesh(
    const Eigen::VectorXd cumsum,
    const double a, 
    const int start, 
    const int stop);


void random_points_on_mesh(
  const int n,
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & X)
{
  X.resize(n,3);
  G generator;
  Real_D distribution(0.0,1.0);

  Eigen::VectorXd double_area;
  igl::doublearea(V, F, double_area);

  Eigen::VectorXd cumsum;
  igl::cumsum(double_area, 1, cumsum);
  // Normalize the results for sampling
  cumsum = cumsum / cumsum(cumsum.rows() - 1);

  for (int i = 0; i < n; i++)
  {
    double a = distribution(generator);
    int idx = triangle_in_mesh(cumsum, a, 0, cumsum.rows());
    Eigen::Vector3d position;
    pt_on_triangle(V.row(F(idx, 0)), 
	           V.row(F(idx, 1)), 
		   V.row(F(idx, 2)), 
		   position);
    X.row(i) = position;
  }
  //for(int i = 0;i<X.rows();i++) X.row(i) = V.row(i%V.rows());
}


////////////////////////////////////////////////////////////////////////////////
// Function bodies
////////////////////////////////////////////////////////////////////////////////
void pt_on_triangle(
    const Eigen::Vector3d v1, 
    const Eigen::Vector3d v2, 
    const Eigen::Vector3d v3, 
    Eigen::Vector3d & position)
{
  G generator;
  Real_D distribution(0.0,1.0);
  double a1 = distribution(generator);
  double a2 = distribution(generator);
  
  if ((a1 + a2) > 1)
  {
    a1 = 1 - a1;
    a2 = 1 - a2;
  }

  position << v1 + \
            a1 * (v2 - v1) + \
            a2 * (v3 - v1); 
}


int triangle_in_mesh(
    const Eigen::VectorXd cumsum,
    const double a, 
    const int start, 
    const int stop)
{
  int left = 0;
  int right = cumsum.rows() - 1;

  while (left <= right)
  {
    int m = (left + right) / 2;
    if (cumsum[m] < a && cumsum[m+1] < a)
    {
      left = m+1;
    }
    else if (cumsum[m] > a && cumsum[m+1] > a)
    {
      right = m-1;
    } else
    {
      return m;
    }
  }
  return 0;
  // if (cumsum(start) > a)
  // {
  //   return start; 
  // }
  // else if (cumsum(start) <= a && cumsum(start + 1) >= a)
  // {
  //   return start + 1; 
  // }

  // int half_size = stop / 2;
  // if (cumsum(half_size) > a)
  // {
  //   return triangle_in_mesh(cumsum, a, start, half_size);
  // } 
  // else
  // {
  //   return triangle_in_mesh(cumsum, a, half_size, stop);
  // }
}
