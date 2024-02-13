#include <assemble_forces.h>
#include <dV_linear_tetrahedron_dq.h>
#include <iostream>

void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::MatrixXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, Eigen::Ref<const Eigen::VectorXd> v0,
                     double C, double D) { 
    f.resize(q.size());
    f.setZero();
    // all tets
    for (int r = 0; r < T.rows(); r ++) {
        Eigen::RowVector4i element = T.row(r);
        Eigen::Vector12d dV;
        dV_linear_tetrahedron_dq(dV, q, V, element, v0(r), C, D);
        // f = - \partial{V}{q} = - \sum \partial{V_r}{q}
        //   = - \sum (E^T {dV})_r
        //  E^T(3 * idx[i] + j, 3i + j) = 1 for i = 0..3, j = 0..2
        //  dV(x) for x = 0..11
        //  (E^T dV)(3 * idx[i] + j) += dV(3i + j) for i = 0..3, j = 0..2
        for (int i = 0; i < 4; i ++) {
            for (int j = 0; j < 3; j ++) {
                f(3*element(i)+j) -= dV(3*i+j);
            }
        }
    }
}