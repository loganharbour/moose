#include "MooseTypes.h"

namespace HCAngularQuadrature
{

/// Provides abscissae and weights of Gauss-Legendre quadrature for x \in [0, 1]
void getLegHRQ(unsigned int order, std::vector<Real> & x, std::vector<Real> & w);

/**
 * Computes the half range Gauss-Legendre-Chebyshev quadrature. Azimuthal order (== cheb_order)
 * and polar order (== leg_order). The cheb_order is the total number of abscissae in phi \in [0,
 * 2 * pi]
 */
void getHalfRangeAQ3D(unsigned int cheb_order,
                      unsigned int leg_order,
                      std::vector<std::pair<Real, Real>> & x,
                      std::vector<Real> & w);
/**
 * Computes the Gauss-Legendre abscissae for a quadrature over polar angle theta for theta \in
 * [-pi / 2, pi / 2]. In 2D it is more advantageous to integrate over angle theta as opposed to mu
 * = cosine(theta). The abscissae are to be interpreted as angles as follows: theta_j =
 * x[j].second, while x[j].first indicates if it theta is positive or negative (mostly unused).
 */
void getHalfRangeAQ2D(unsigned int leg_order,
                      std::vector<std::pair<Real, Real>> & x,
                      std::vector<Real> & w);

/// Returns a vector orthogonal to v
Point getOrthonormalVector(const Point & v, unsigned int dim);

/**
 * The angular quadrature can be used to create angular directions w.r.t. to
 * a reference vector, e.g. the unit vector in z -- ez. For the view factor algorithm
 * we need the direction to be rotated to be in the coordinate system of the local
 * normal. This method computes the rotation matrix to accomplish that.
 */
DenseMatrix<Real> aqRoationMatrix(const Point & normal, const unsigned int dim);

Point getAngularDirection(const unsigned int l,
                          const DenseMatrix<Real> & rotation_matrix,
                          const std::vector<std::pair<Real, Real>> & angles,
                          const unsigned int dim);

}
