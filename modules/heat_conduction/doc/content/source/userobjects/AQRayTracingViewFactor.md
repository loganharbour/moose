# AQRayTracingViewFactor

## Description

`AQRayTracingViewFactor` uses ray tracing for computing view factors in two-dimensional and three-dimensional cavities. It does not impose
any restrictions on the geometry of the cavity; in particular it allows non-planar surfaces in radiative exchange, non-convex cavities, and obstruction
of view between two surfaces by a third surface.

The central idea for computing view factors using ray tracing is to transform the integral over the target area into an integral over
angular direction (i.e., an integral over the field of view of any infinitesimal element on the starting surface). An angular quadrature
is used to numerically approximate the angular integral. To this end, the ray tracing module is used to follow rays along the directions
of the angular quadrature and determine which surface they intersect first. The ray is terminated on that surface and the contribution
to the view factor between the surface of origin and this surface is incremented.

The view factor from the surface indexed by $i$ to the surface indexed by $j$, $F_{i,j}$, is defined as the following double integral:

\begin{equation}\label{eq:view_factor}
  F_{i,j} = \frac{1}{A_i \pi} \int_{A_i} \int_{A_j}  \frac{\cos \beta_i \cos \beta_j}{r^2}  dA_i dA_j,
\end{equation}

where $A_i$ is the area of surface $i$, $r$ is the distance between the centroids of the infinitesimal areas $dA_i$ and $dA_j$,
and $\beta_i$ is the angle between the normal at $dA_i$ and the line connecting $dA_i$ and $dA_j$. The solid angle $d \hat{\Omega}$
that is subtended by the infinitesimal element $d A_j$ when observed from the centroid of $dA_i$ is:

\begin{equation}\label{eq:trafo}
  dA_j = \frac{r^2}{\cos \beta_j} d \hat{\Omega}.
\end{equation}

Using [eq:trafo], [eq:view_factor] is rewritten as:

\begin{equation}\label{eq:trafo_view_factor}
 F_{i,j} = \frac{1}{A_i \pi} \int_{A_i} \int_{\hat{\Omega}_j} \hat{\Omega}\cdot \vec{n}  ~dA_i d\hat{\Omega},
\end{equation}

where $\hat{\Omega}_j$ is the solid angle that surface $j$ subtends without taking into considerations obstructions, $\vec{n}$ is the normal at the centroid of $dA_i$ and $\cos \beta_i=\hat{\Omega}\cdot \vec{n}$.

Before discretizing [eq:trafo_view_factor], we rewrite it slightly by introducing the function $H_j\left(\vec{r},\hat{\Omega}\right)$ which is $1$ if
a ray starting from location $\vec{r}$ into direction $\hat{\Omega}$ makes its first intersection with surface $j$. Then [eq:trafo_view_factor] becomes

\begin{equation}\label{eq:trafo_view_factor_2}
 F_{i,j} = \frac{1}{A_i \pi} \int_{A_i} \int_{2 \pi^+} \hat{\Omega}\cdot \vec{n} H_j\left(\vec{r},\hat{\Omega}\right) ~dA_i d\hat{\Omega},
\end{equation}

where the angular integral extends over half of the angular range directed into the cavity denoted by $2 \pi^+$.
The integral in [eq:trafo_view_factor_2] is discretized using a spatial quadrature over the area $i$ and an angular quadrature over $2 \pi^+$:

\begin{equation}\label{eq:trafo_trafo}
 F_{i,j} = \frac{1}{A_i \pi} \sum\limits_{l=1}^L  \sum\limits_{k=1}^K w_l \omega_k  \left(\hat{\Omega}_k\cdot \vec{n}_l\right) H_j\left(\vec{r}_l,\hat{\Omega}_k\right) ,
\end{equation}

where $l$ enumerates the spatial quadrature points, $k$ enumerates the angular directions, $\vec{r}_l$ is the location of spatial quadrature point $l$, $\vec{n}_l$ is the
normal at spatial quadrature point $l$, and $\hat{\Omega}_k$ is the k-th angular direction.

`AQRayTracingViewFactor` uses Gaussian product quadratures for integrals in space, and a half-range Gauss-Chebyshev quadrature in angle. The Gauss-Chebyshev is
identical to the standard full range Gauss-Chebyshev quadrature except that the polar integral is taken only from $0$ to $\pi / 2$ instead of $-\pi/2$ to $\pi / 2$.

Symmetry surfaces require special treatment. The difference between surfaces that participate in the radiative transfer and symmetry surfaces is that
view factors are not computed for symmetry surfaces (i.e. if either surface $i$ or $j$ in $F_{i,j}$ is a symmetry surface, the view factor is not computed).
Rays neither start nor end on symmetry surfaces. Instead, rays are specularly reflected off of symmetry surfaces. This is facilitated by the ray tracing
module by using the `ReflectRayBC`.  

## Example Input syntax

!listing modules/heat_conduction/test/tests/view_factors/analytical_with_obstruction_aq.i
block=UserObjects

!syntax parameters /UserObjects/AQRayTracingViewFactor

!syntax inputs /UserObjects/AQRayTracingViewFactor

!syntax children /UserObjects/AQRayTracingViewFactor
