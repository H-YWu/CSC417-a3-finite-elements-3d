#include <build_skinning_matrix.h>
#include <phi_linear_tetrahedron.h>
#include <vector>

// Utils: check if the point X is at the same side with another vertex to the corresponding surface 
bool same_side_tet(Eigen::Ref<const Eigen::Vector3d> V0,
                   Eigen::Ref<const Eigen::Vector3d> V1,
                   Eigen::Ref<const Eigen::Vector3d> V2,
                   Eigen::Ref<const Eigen::Vector3d> V3,
                   Eigen::Ref<const Eigen::Vector3d> X) {
    Eigen::Vector3d V01 = V1 - V0;
    Eigen::Vector3d normal = V01.cross(V2-V0);
    double dotV3 = normal.dot(V3-V0);
    double dotX = normal.dot(X-V0);
    if (dotV3 <= 0.0 && dotX <= 0.0) return true;
    if (dotV3 >= 0.0 && dotX >= 0.0) return true;
    return false;
}

// Utils: check if the point X is inside the tet
bool point_in_tet(Eigen::Ref<const Eigen::Vector3d> V0,
                  Eigen::Ref<const Eigen::Vector3d> V1,
                  Eigen::Ref<const Eigen::Vector3d> V2,
                  Eigen::Ref<const Eigen::Vector3d> V3,
                  Eigen::Ref<const Eigen::Vector3d> X) {
    return same_side_tet(V0, V1, V2, V3, X) &&
           same_side_tet(V1, V2, V3, V0, X) &&
           same_side_tet(V2, V3, V0, V1, X) &&
           same_side_tet(V3, V0, V1, V2, X);
}

void build_skinning_matrix(Eigen::SparseMatrixd &N, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, 
                                                   Eigen::Ref<const Eigen::MatrixXd> V_skin) {
    N.resize(V_skin.rows(), V.rows());
    std::vector<Eigen::Triplet<double>> triples;
    // X_{surface} = N * V
    //  lx3 = lxn x nx3
    //  X_{surface}.row(i) = N.row(i) * V
    for (int i = 0; i < V_skin.rows(); i ++) {
        Eigen::Vector3d X = V_skin.row(i);
        // Check all tets to find which one contains the X
        for (int j = 0; j < T.rows(); j ++) {
            Eigen::RowVector4i element = T.row(j);
            Eigen::Vector3d V0 = V.row(element(0));
            Eigen::Vector3d V1 = V.row(element(1));
            Eigen::Vector3d V2 = V.row(element(2));
            Eigen::Vector3d V3 = V.row(element(3));
            if (point_in_tet(V0, V1, V2, V3, X)) {
                Eigen::Vector4d phi;
                phi_linear_tetrahedron(phi, V, element, X);
                //  N[i, element(k)] = phi(Xi)[k] for k=0..3
                for (int k = 0; k < 4; k ++) {
                    triples.push_back(Eigen::Triplet<double>(i, element(k), phi(k)));
                }
                break;
            }
        }
    }
    N.setFromTriplets(triples.begin(), triples.end());
}