#include <matrix_BE_tetrahedron.h>

void matrix_B_tetrahedron(Eigen::SparseMatrixd &B, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element) {
    B.resize(9, 12);
    // the 3x3 delta X matrix of a single tetrahedron 
    Eigen::Matrix3d T;
    std::vector<Eigen::Vector3d> X(4);
    for (int i = 0; i < 4; i ++) {
        X[i] = V.row(element(0, i));
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
    // B(j + 3k, 3i + k) = D(i, j) for i = 0..3, j = 0..2, k = 0..2
    std::vector<Eigen::Triplet<double>> triples;
    for (int i = 0; i < 4; i ++) {
        for (int j = 0; j < 3; j ++) {
            for (int k = 0; k < 3; k ++) {
                triples.push_back(Eigen::Triplet<double>(j+k*3, i*3+k, D(i,j)));
            }
        }
    }
    B.setFromTriplets(triples.begin(), triples.end());
}