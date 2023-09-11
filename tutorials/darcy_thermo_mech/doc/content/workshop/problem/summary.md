# Tutorial Steps

!style halign=center
[#step01](problem/step01.md)\\
[#step02](problem/step02.md)\\
[#step03](problem/step03.md)\\
[#step04](problem/step04.md)\\
[#step05](problem/step05.md)\\
[#step06](problem/step06.md)\\
[#step07](problem/step07.md)\\
[#step08](problem/step08.md)\\
[#step09](problem/step09.md)\\
[#step10](problem/step10.md)\\
[#step11](problem/step11.md)

!---

## [#step01](problem/step01.md)

The first step is to solve a simple "Diffusion" problem, which requires no code. This step
will introduce the basic system of MOOSE.

!equation
-\nabla\cdot\nabla p = 0

!---

## [#step02](problem/step02.md)

In order to implement the Darcy pressure equation, a `Kernel` object is needed to represent:

!equation
-\nabla\cdot\frac{\mathbf{K}}{\mu}\nabla p = 0

!---

## [#step03](problem/step03.md)

Instead of passing constant parameters to the pressure diffusion `Kernel` object, the Material
system can be used to supply the values. This allows for properties that vary in space and time
as well as be coupled to variables in the simulation.

!---

## [#step04](problem/step04.md)

The velocity is computed from the pressure based on Darcy's law as:

!equation
\vec{u} = -\frac{\mathbf{K}}{\mu} \nabla p

This velocity can be computed using the Auxiliary system.

!---

## [#step05](problem/step05.md)

Solve the transient heat equation using the "heat conduction" module.

!equation
C \frac{\partial T}{\partial t} - \nabla \cdot k \nabla T = 0

!---

## [#step06](problem/step06.md)

Solve the pressure and temperature in a coupled system of equations by adding the advection term
to the heat equation.

!equation
-\nabla \cdot \frac{\mathbf{K}}{\mu} \nabla p  = 0
\\
C\left( \frac{\partial T}{\partial t} + \epsilon \vec{u}\cdot\nabla T \right) - \nabla \cdot k \nabla T = 0

!---

## [#step07](problem/step07.md)

In the transient simulation, a "traveling wave" profile moves through the porous medium. Instead of
using a uniform mesh to resolve the wave profile, we can dynamically adapt the mesh to the solution.

!---

## [#step08](problem/step08.md)

`Postprocessor` and `VectorPostprocessor` objects can be used to compute aggregate value(s) for a
simulation, such as the average temperature on the boundary or the temperatures along a line
within the solution domain.

!---

## [#step09](problem/step09.md)

Thermal expansion of the porous media can be added to the coupled set of equations
using the "tensor mechanics" module, without adding additional code.

!---

## [#step10](problem/step10.md)

MOOSE is capable of running multiple applications together and transfer data between the various
applications.

This problem replaces the thermal conductivity calculated by the Material with a value computed by
another application that runs a phase-based micro-structure simulation.

!---

## [#step11](problem/step11.md)

MOOSE includes a system to create custom input syntax for common tasks, in this step the syntax
for the two equations and velocity auxiliary calculation are simplified for end-users.
