#include <mass_matrix_mesh.h>
#include <mass_matrix_linear_tetrahedron.h>
#include <iostream>

void mass_matrix_mesh(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::MatrixXi> T, double density, Eigen::Ref<const Eigen::VectorXd> v0) {
    M.resize(qdot.size(), qdot.size());
    // all tets
    std::vector<Eigen::Triplet<double>> triples;
    for (int r = 0; r < T.rows(); r ++) {
        Eigen::RowVector4i element = T.row(r);
        Eigen::Matrix1212d m;
        mass_matrix_linear_tetrahedron(m, qdot, element, density, v0(r));
        // E(3i + j, 3 * idx[i] + j)   = 1 for i = 0..3, j = 0..2
        // E^T(3 * idx[i] + j, 3i + j) = 1 for i = 0..3, j = 0..2
        // m(x, y) for x,y = 0..11
        //  (m E)(x, 3 * idx[i] + j) += m(x, 3i + j) for x = 0..11, i = 0..3, j = 0..2
        // M: (E^T m E) (3 * idx[i] + j, 3 * idx[i] + j) += m(3i + j, 3i + j) for i = 0..3, j = 0..2
        for (int i = 0; i < 4; i ++) {
            for (int j = 0; j < 3; j ++) {
                triples.push_back(Eigen::Triplet<double>(3*element(i)+j, 3*element(i)+j, m(3*i+j,3*i+j)));
            }
        }

    }
    M.setFromTriplets(triples.begin(), triples.end());
}