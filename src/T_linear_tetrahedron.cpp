#include <T_linear_tetrahedron.h>
#include <mass_matrix_linear_tetrahedron.h>
#include <iostream>

void T_linear_tetrahedron(double &T, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::RowVectorXi> element, double density, double volume) {
    Eigen::Matrix1212d M;
    mass_matrix_linear_tetrahedron(M, qdot, element, density, volume);
    // T = 1/2 (E_j qdot)^T M_j (E_j qdot)
    Eigen::VectorXd element_qdot;
    for (int i = 0; i < 4; i ++) {
        element_qdot.segment(i * 3, 3) = qdot.segment(element(i) * 3, 3);
    }
    T = 0.5 * element_qdot.transpose() * M * element_qdot;
}