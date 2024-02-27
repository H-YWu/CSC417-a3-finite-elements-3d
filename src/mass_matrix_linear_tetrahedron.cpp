#include <mass_matrix_linear_tetrahedron.h>

void mass_matrix_linear_tetrahedron(Eigen::Matrix1212d &M, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::RowVectorXi> element, double density, double volume) {
    Eigen::MatrixXd Me(4, 4);
    // generated using MATLAB symbolic toolkit
    Me(0,0) = 1.0/1.0E+1;
    Me(0,1) = 1.0/2.0E+1;
    Me(0,2) = 1.0/2.0E+1;
    Me(0,3) = 1.0/2.0E+1;
    Me(1,0) = 1.0/2.0E+1;
    Me(1,1) = 1.0/1.0E+1;
    Me(1,2) = 1.0/2.0E+1;
    Me(1,3) = 1.0/2.0E+1;
    Me(2,0) = 1.0/2.0E+1;
    Me(2,1) = 1.0/2.0E+1;
    Me(2,2) = 1.0/1.0E+1;
    Me(2,3) = 1.0/2.0E+1;
    Me(3,0) = 1.0/2.0E+1;
    Me(3,1) = 1.0/2.0E+1;
    Me(3,2) = 1.0/2.0E+1;
    Me(3,3) = 1.0/1.0E+1;
    M.setZero();
    for (int i = 0; i < 4; i ++) {
        for (int j = 0; j < 4; j ++) {
            double e = density * volume * Me(i, j);
            for (int k = 0; k < 3; k ++) {
                M(i * 3 + k, j * 3 + k) = e; 
            }
        }
    }  
}