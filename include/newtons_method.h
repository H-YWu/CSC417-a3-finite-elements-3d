#include <Eigen/Dense>
#include <EigenTypes.h>

//Input:
//  x0 - initial point for newtons search
//  f(x) - function that evaluates and returns the cost function at x
//  g(dx, x) - function that evaluates and returns the gradient of the cost function in dx
//  H(dH, x) - function that evaluates and returns the Hessian in dH (as a sparse matrix).
//  max steps - the maximum newton iterations to take
//  tmp_g and tmp_H are scratch space to store gradients and hessians
//Output: 
//  x0 - update x0 to new value
template<typename Objective, typename Jacobian, typename Hessian>
double newtons_method(Eigen::VectorXd &x0, Objective &f, Jacobian &g, Hessian &H, unsigned int maxSteps, Eigen::VectorXd &tmp_g, Eigen::SparseMatrixd &tmp_H) {
	// Newton's Method with backtracking line search
	/**
	* i = 0; x0 = initial guess;
	* while i < max steps do:
	*    // check convergence
	*    if || g^i || < tol: break; 
	*    // search direction
	*    d := solve linear system H^i d = - g^i;
	*    // chosse alpha: backtracking line search
	*    alpha = alpha_max; find = false; p < 1.0; c = choose;
	*    while not find do:
	*       if (cost(x0<i> + alpha * d) <= cost(x0<i>) + c * d^T g^i) ||
	*          alpha < tol: find = true;
	*       else: alpha = p * alpha;
	*    end
	*    // update
	*    x0<i+1> = x0<i> + alpha * d; i = i + 1;
	* end
	*/

	// constants
	double tol = 1e-4, tol2 = 1e-8, max_alpha = 1.0, p = 0.5, c = 1e-8;

	Eigen::VectorXd d;
	for (int i = 0; i < maxSteps; i ++) {
	    // check convergence
		g(tmp_g, x0);
		if (tmp_g.dot(tmp_g) < tol2) break;
		H(tmp_H, x0);
		// search direction
		//	solve H^i d = -g^i (dH d = -dx) 
		Eigen::SimplicialLDLT<Eigen::SparseMatrixd> solver;
		solver.compute(tmp_H);
		if(solver.info() != Eigen::Success) {
			std::cerr << "ERROR: H decomposition failed\n";
			return -1.0;
		}
		d = solver.solve(-tmp_g);
		if(solver.info() != Eigen::Success) {
			std::cerr << "ERROR: Solving Hd = -g failed\n";
			return -1.0;
		}
	    // chosse alpha: backtracking line search
		double alpha = max_alpha, threshold = f(x0) + c * d.dot(tmp_g);
		while (true) {
			if ((f(x0 + alpha * d) <= threshold) ||
				alpha < tol) break;
			alpha *= p;
		}
		// update
		x0 = x0 + alpha * d;
	}

	return 0.0;
}
