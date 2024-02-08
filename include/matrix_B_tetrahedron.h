#include <Eigen/Dense>
#include <EigenTypes.h>

//Input:
// V - vertex matrix for the mesh
// element - vertex indices of the element
//Output:
// B - the 9x12 B auxilary matrix of a single tetrahedron to calculate stacked F
void matrix_B_tetrahedron(Eigen::SparseMatrixd &B, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element);