#include <d2V_spring_particle_particle_dq2.h>

void d2V_spring_particle_particle_dq2(Eigen::Ref<Eigen::Matrix66d> H, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness) {
    // q = [q0; q1]
    // dx = Bq (q0 -> q1)
    //  B = [-I, I]
    // f = - 1/2 k (\sqrt{dx^T dx} - l0) \frac{dx^T}{\sqrt{dx^T dx}} B
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    Eigen::Matrix66d BTB;
    BTB << I, -I,
            -I, I;
    Eigen::Vector3d dx = q1 - q0;
    // transpose
    Eigen::Vector6d BTdx;
    BTdx << -dx, dx;
    double dxdot = dx.dot(dx);
    // (+) H = \partial^2{V}{q^2} = -\partial{f}{q}
    //       = -\partial{f}{dx} \diff{dx}{q}
    //       = 1/2 k (\frac{\sqrt{dx^T dx} - l0}{\sqrt{dx^T dx}} B^T + \frac{l0}{2 (dx^T dx)^{3/2}} B^T dx dx^T) B
    //       = 1/2 k (\frac{\sqrt{dx^T dx} - l0}{\sqrt{dx^T dx}} B^T B + \frac{l0}{2 (dx^T dx)^{3/2}} B^T dx dx^T B)
    //       = 1/2 k (\frac{\sqrt{dx^T dx} - l0}{\sqrt{dx^T dx}} B^T B + \frac{l0}{2 (dx^T dx)^{3/2}} B^T dx (B^T dx)^T)
    H = 0.5 * stiffness * 
        ( ((sqrt(dxdot) - l0) * pow(dxdot, -0.5)) * BTB +
          (0.5 * l0 * pow(dxdot, -1.5)) * BTdx * BTdx.transpose());
}