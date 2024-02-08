#include <dphi_linear_tetrahedron_dX.h>
#include <phi_linear_tetrahedron.h>
#include <iostream>
void dphi_linear_tetrahedron_dX(Eigen::Matrix43d &dphi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
    // lambda = T^{-1} * b
    // \partial{lambda}{X} = \partial{T^{-1} * b}{X}
    //                     = \partial{T^{-1} * b}{b} \diff{b}{X}
    //                     = \partial{T^{-1} * b}{b} \diff{X - X0}{X}
    //                     = T^{-1} I
    //                     = T^{-1}
    // phi = [1 - lambda[+012]; lambda]
    // \partial{phi}{X} = [-\partial{lambda}{X}[+012]; \partial{lambda}{X}]
    std::vector<Eigen::Vector3d> XX(4);
    for (int i = 0; i < 4; i ++) {
        XX[i] = V.row(element(0, i));
    }
    Eigen::MatrixXd T(3, 3);
    for (int i = 0; i < 3; i ++) {
        T.col(i) = XX[i+1] - XX[0];
    }
    Eigen::MatrixXd Tinv = T.inverse();
    dphi.row(0) << 0., 0., 0.;
    for (int i = 0; i < 3; i ++) {
        dphi.row(i+1) = Tinv.row(i);
        dphi.row(0) -= Tinv.row(i);
    }
}