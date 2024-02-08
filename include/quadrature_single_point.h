#include <Eigen/Dense>
#include <EigenTypes.h>

//Input:
//  q - generalized coordinates of the FEM system
//  element - vertex indices for the tetrahedron
// volume - volume of the tetrahedron
// integrand(out, q, X) - function to be integrated, returns value in out.
//Output:
//  integrated - the value of the integrated function
template<typename Ret, typename Integrand_Func>
inline void quadrature_single_point(Ret &&integrated, Eigen::Ref<const Eigen::VectorXd> q, 
                                               Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                                               Integrand_Func integrand) {
    // Take the centroid as the single quadrature point
    //  However it's actually useless in the linear case
    Eigen::Vector3d X;
    X << 0.0, 0.0, 0.0;
    for (int i = 0; i < 4; i ++) {
        X += 0.25 * q.segment(element(i) * 3, 3);
    }
    integrand(integrated, q, element, X);
    integrated = volume * integrated;
}

