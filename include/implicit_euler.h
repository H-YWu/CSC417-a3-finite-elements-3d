#include <Eigen/Dense>
#include <EigenTypes.h>
#include <newtons_method.h>

//Input:
//  q - generalized coordinates for the FEM system
//  qdot - generalized velocity for the FEM system
//  dt - the time step in seconds
//  mass - the mass matrix
//  energy(q, qdot) -  a function that computes the energy of the FEM system. This takes q and qdot as parameters, returns the energy value.
//  force(f, q, qdot) - a function that computes the force acting on the FEM system. This takes q and qdot as parameters, returns the force in f.
//  stiffness(K, q, qdot) - a function that computes the stiffness (negative second derivative of the potential energy). This takes q and qdot as parameters, returns the stiffness matrix in K.  
//  tmp_qdot - scratch space for storing velocities 
//  tmp_force - scratch space to collect forces
//  tmp_stiffness - scratch space to collect stiffness matrix
//Output:
//  q - set q to the updated generalized coordinate using linearly implicit time integration
//  qdot - set qdot to the updated generalized velocity using linearly implicit time integration
template<typename ENERGY, typename FORCE, typename STIFFNESS> 
inline void implicit_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, 
                            const Eigen::SparseMatrixd &mass,  ENERGY &energy, FORCE &force, STIFFNESS &stiffness, 
                            Eigen::VectorXd &tmp_qdot, Eigen::VectorXd &tmp_force, Eigen::SparseMatrixd &tmp_stiffness) {
    tmp_qdot = qdot;

    // gradient of the cost energy
    auto g = [&](Eigen::VectorXd &dx, Eigen::Ref<const Eigen::VectorXd> x) { 
        // dx := M qdot^{t+1}_{i} − M qdot^{t} + dt * \partial{V(q^{t} + dt * qdot^{t+1}_{i})}{q}
        // x := qdot&{t+1}_{i}
        force(dx, q + dt * x, x);
        dx *= -dt;
        dx += mass * (x - qdot);
    };

    // hessian of the cost energy
    auto H = [&](Eigen::SparseMatrixd &dH, Eigen::Ref<const Eigen::VectorXd> x) {
        // dH := M + dt^2 * \partial^2{V(q^{t} + dt * qdot^{t+1}_{i})}{q}
        // x := qdot&{t+1}_{i}
        stiffness(dH, q + dt * x, x);
        dH *= (dt*dt);
        dH += mass;
    };

    // To ensure reasonable performance
    int maxIters = 5;

    newtons_method(tmp_qdot, energy, g, H, maxIters, tmp_force, tmp_stiffness);
    qdot = tmp_qdot;
    q = q + dt * qdot;
}
