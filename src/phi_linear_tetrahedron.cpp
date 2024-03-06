#include <phi_linear_tetrahedron.h>

void phi_linear_tetrahedron(Eigen::Vector4d &phi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
    std::vector<Eigen::Vector3d> XX(4);
    for (int i = 0; i < 4; i ++) {
        XX[i] = V.row(element(i));
    }
    Eigen::MatrixXd T(3, 3);
    for (int i = 0; i < 3; i ++) {
        T.col(i) = XX[i+1] - XX[0];
    }
    Eigen::Vector3d b = X - XX[0];
    Eigen::Vector3d lamb = T.inverse() * b;
    phi.coeffRef(0) = 1.0;
    for (int i = 0; i < 3; i ++) {
        phi.coeffRef(i+1) = lamb.coeff(i);
        phi.coeffRef(0) -= lamb.coeff(i);
    }
    Eigen::Vector3d test;
    test << 0.0, 0.0, 0.0;
    for (int i = 0; i < 4; i ++) {
        test += phi(i) * XX[i];
    }
}