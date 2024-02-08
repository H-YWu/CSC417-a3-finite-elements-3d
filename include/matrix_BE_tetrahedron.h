#include <Eigen/Dense>
#include <EigenTypes.h>

//Input:
// q - generalized coordinates of the FEM system
// V - vertex matrix for the mesh
// element - vertex indices of the element
//Output:
// BE - the 9x3n B_j * E_j auxilary matrix of a single tetrahedron, which is \partial{F}{q}
void matrix_BE_tetrahedron(Eigen::SparseMatrixd &BE,  Eigen::Ref<const Eigen::VectorXd> q, 
                        Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element);