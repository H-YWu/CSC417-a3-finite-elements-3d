#include <dV_linear_tetrahedron_dq.h>

#include <dphi_linear_tetrahedron_dX.h>
#include <dpsi_neo_hookean_dF.h>
#include <quadrature_single_point.h>
#include <iostream>
#include <matrix_B_tetrahedron.h>
#include <deformation_gradient.h>

void dV_linear_tetrahedron_dq(Eigen::Vector12d &dV, Eigen::Ref<const Eigen::VectorXd> q, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                          double C, double D) {

   auto neohookean_linear_tet = [&](Eigen::Vector12d &dV, Eigen::Ref<const Eigen::VectorXd>q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
        Eigen::Matrix3d F;
        deformation_gradient(F, q, V, element);
        Eigen::Vector9d dpsi;
        dpsi_neo_hookean_dF(dpsi, F, C, D); 
        Eigen::SparseMatrixd B;
        matrix_B_tetrahedron(B, V, element);
        Eigen::SparseMatrixd BT = B.transpose();
        dV = BT * dpsi;
    };

    quadrature_single_point(dV, q, element, volume, neohookean_linear_tet);  
}