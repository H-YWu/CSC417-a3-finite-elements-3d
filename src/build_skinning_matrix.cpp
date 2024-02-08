#include <build_skinning_matrix.h>
#include <phi_linear_tetrahedron.h>
#include <vector>
#include <iostream>

void build_skinning_matrix(Eigen::SparseMatrixd &N, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, 
                                                   Eigen::Ref<const Eigen::MatrixXd> V_skin) {
    N.resize(V_skin.rows(), V.rows());
    std::vector<Eigen::Triplet<double>> triples;
    for (int i = 0; i < V_skin.rows(); i ++) {
        Eigen::Vector3d X = V_skin.row(i);
        for (int j = 0; j < T.rows(); j ++) {
            Eigen::RowVector4i element = T.row(j);
            Eigen::Vector4d phi;
            phi_linear_tetrahedron(phi, V, element, X);
            bool find = true;
            for (int k = 0; k < 4; k ++) {
                if (phi(k) < 0.0) {
                    find = false;
                    break;
                }
            }
            if (find) {
                for (int k = 0; k < 4; k ++) {
                    triples.push_back(Eigen::Triplet<double>(i, element(k), phi(k)));
                }
                break;
            }
        }
    }
    N.setFromTriplets(triples.begin(), triples.end());
}