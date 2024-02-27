#include <deformation_gradient.h>

void deformation_gradient(Eigen::Matrix3d &F, Eigen::Ref<const Eigen::VectorXd> q, 
                        Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element){
    // the 3x3 delta X matrix of a single tetrahedron 
    Eigen::Matrix3d T;
    std::vector<Eigen::Vector3d> X(4);
    Eigen::Matrix34d x;
    for (int i = 0; i < 4; i ++) {
        X[i] = V.row(element(i));
        x.col(i) = q.segment(element(i)*3, 3);
    }
    for (int i = 0; i < 3; i ++) {
        T.col(i) = X[i+1] - X[0];
    }
    Eigen::RowVector3d neg1T;
    neg1T << -1., -1., -1.;
    Eigen::Matrix3d Tinv = T.inverse();
    // auxilary matrix D = [-\vec{1}^T T^{-1}; T^{-1}]
    Eigen::Matrix43d D;
    D.row(0) = neg1T * Tinv;
    D.block(1, 0, 3, 3) = Tinv;
    // F = [x0, x1, x2, x3] * [-\vec{1}^T T^{-1}; T^{-1}]
    F = x * D;
}