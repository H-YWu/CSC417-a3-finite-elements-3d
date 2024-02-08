#include <psi_neo_hookean.h>
#include <dphi_linear_tetrahedron_dX.h>
#include <iostream>
void psi_neo_hookean(double &psi, 
                     Eigen::Ref<const Eigen::Matrix3d> F,
                     double C, double D) {
    double J = F.determinant();
    Eigen::MatrixXd FTF = F.transpose() * F;
    psi = C * (pow(J, -2/3) * FTF.trace() - 3.) + D * pow(J - 1., 2);
}