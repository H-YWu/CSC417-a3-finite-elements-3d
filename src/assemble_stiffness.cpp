#include <assemble_stiffness.h>
#include <d2V_linear_tetrahedron_dq2.h>

void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, Eigen::Ref<const Eigen::VectorXd> v0, 
                     double C, double D) { 
    K.resize(qdot.size(), qdot.size());
    // all tets
    std::vector<Eigen::Triplet<double>> triples;
    for (int r = 0; r < T.rows(); r ++) {
        Eigen::RowVector4i element = T.row(r);
        Eigen::Matrix1212d H;
        d2V_linear_tetrahedron_dq2(H, q, V, element, v0(r), C, D); 
        // K = - \partial^2{V}{q} = - \sum \partial^2{V_r}{q}
        //   = - \sum (E^T {H} E)_r
        //  E(3i + j, 3 * idx[i] + j)   = 1 for i = 0..3, j = 0..2
        //  E^T(3 * idx[i] + j, 3i + j) = 1 for i = 0..3, j = 0..2
        //  H(x, y) for x,y = 0..11
        //  (H E)(x, 3 * idx[i] + j) -= H(x, 3i + j) for x = 0..11, i = 0..3, j = 0..2
        //  (E^T H E)(3 * idx[I] + J, 3 * idx[i] + j) -= H(3I + J, 3i + j) for I = 0..3, J = 0..2, i = 0..3, j = 0..2
        for (int I = 0; I < 4; I ++) {
            for (int J = 0; J < 3; J ++) {
                for (int i = 0; i < 4; i ++) {
                    for (int j = 0; j < 3; j ++) {
                        triples.push_back(Eigen::Triplet<double>(3*element(I)+J, 3*element(i)+j, -H(3*I+J,3*i+j)));
                    }
                }
            }
        }
    }
    K.setFromTriplets(triples.begin(), triples.end());
};
