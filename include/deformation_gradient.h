#include <Eigen/Dense>
#include <EigenTypes.h>

//Input:
// q - generalized coordinates of FEM system
// V - vertex matrix for the mesh
// element - vertex indices of the element
//Output:
// F - the 3x3 deformation gradient of a single tetrahedron
void deformation_gradient(Eigen::Matrix3d &F, Eigen::Ref<const Eigen::VectorXd> q, 
                        Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element);