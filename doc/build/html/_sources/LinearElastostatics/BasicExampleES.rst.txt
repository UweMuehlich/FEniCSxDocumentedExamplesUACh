.. _Target_MinimalExamplePES:

Minimal example for plane elastostatics
---------------------------------------

This is a minimal example for plane elasticity. It will be used as point of departure for the
remaining examples in this part.

We consider a plane problem of linear elasticity. Considering a square with given dimension  :math:`h`, i.e.,  a domain
:math:`\Omega = (0,h) \times (0,h)`,
the governing field equations in  :math:`\Omega` are

.. math::
   \begin{eqnarray}
     \nabla \cdot \sigma + \mathbf{f} &=& 0 \\
                  e(\mathbf{u}) &=& \frac{1}{2} \left[ \nabla \mathbf{u} + (\nabla \mathbf{u})^\mathrm{T} \right] \\
                  \sigma &=& C : e(\mathbf{u})
   \end{eqnarray}

with stress tensor :math:`\sigma`, strain tensor  :math:`e`, displacement vector :math:`\mathbf{u}`, elasticity
tensor :math:`C`, and volume force vector :math:`\mathbf{f}`.

In view of problems to be addressed later on, the square consists of three horizontal stripes. At the bottom
:math:`\partial \Omega_\mathrm{B}`, the square is completely fixed, whereas at the
left boundary (:math:`\partial \Omega_\mathrm{L}`), only the horizontal displacement is blocked. In addition,
a vertical displacement :math:`\bar{u}_2` is prescribed at the top (:math:`\partial \Omega_\mathrm{T}`).
Furthermore, a volume force can be applied.


The corresponding variational problem reads

.. math::
   \begin{eqnarray}
     \int \limits_\Omega  \sigma : e(\mathbf{v}) \, \mathrm{d}^2 x -  \int \limits_\Omega  {f} \cdot \mathbf{v}
     \, \mathrm{d}^2 x   &=& 0 \\
                  \mathbf{u} &=& \mathbf{0} \quad \mathrm{~at~} \partial \Omega_{\mathrm{B}} \\
                  u_1 &=& 0 \quad \mathrm{~at~} \partial \Omega_{\mathrm{L}} \\
                  u_2 &=& \bar{u}_2 \quad \mathrm{at~} \partial \Omega_{\mathrm{T}} \\
   \end{eqnarray}

with test function :math:`\mathbf{v} \in V` and :math:`\mathbf{v}=\mathbf{0}` at
:math:`\partial  \Omega_\mathrm{B}`,  :math:`v_1 = 0` at :math:`\partial  \Omega_\mathrm{L}`
as well as :math:`v_2 = 0` at :math:`\partial  \Omega_\mathrm{B}`, but otherwise arbitrary.
Since, :math:`\mathbf{u}` and :math:`\mathbf{v}` are vector functions, :math:`V` is a vector function space.

.. note::
   You can download the complete python script
   :download:`homPlaneElast.py <homPlaneElast.py>`

In the following, the main parts of the script are discussed in more detail.
First, the necessary modules are imported.  Different ways to import module or just specific functions
of modules exist. Here, `explicit imports
<https://betterprogramming.pub/namespacing-with-python-79574d125564>`_  are preferred. The DOLFINx version
is printed, because there can be changes regarding the implementation depending on the version. At the
moment DOLFINx version: 0.8.0 is used.

.. literalinclude:: homPlaneElast.py
   :lines: 1-17

Next, some general variables are defined, among them the names for the files used for storing
the results.

.. literalinclude:: homPlaneElast.py
   :lines: 18-30

A mesh of quadrilateral elements defined on a rectangular domain is generated using a built-in function.
Afterwards, a  compatible vector function space :math:`V` is defined. The vector character is explicit due to (gdim, ),
since gdim=2.

.. literalinclude:: homPlaneElast.py
   :lines: 31-42


Implementing homogeneous and inhomogeneous Dirichlet boundary conditions involves a number of steps:

* Define the values for the components of the displacement vector to be prescribed.
* Define functions to identify the corresponding boundary nodes (vertices) by means of a criterion.
* Pass through the boundary facets (here edges) checking if the vertices of a facet fulfill the criterion.
* Identify the corresponding degrees of freedom (dof's). **If only specific displacement components should
  be prescribed, the corresponding subspaces must be specified**.
* Assign the prescribed values to the dof's, specifying again the corresponding subspaces if necessary.
* Join all Dirichlet boundary conditions in a list.

.. literalinclude:: homPlaneElast.py
   :lines: 43-75

Constituve relations are defined for plane stress and plane strain, working eventually with Lame's
constants based on Young's modulus and Poisson's number.

.. literalinclude:: homPlaneElast.py
   :lines: 76-96

The definition of the variational problem is rather straight forward, since  the syntax is almost
identical with the FEniCs syntax.

.. literalinclude:: homPlaneElast.py
   :lines: 97-109

After solving the problem, the results can be written into one or more files to be processed, for
instantce with ParaView. Displacement results are given at the nodes of the mesh. Therefore, they
can be written directly into the result file(s). Stresses on the other hand are computed at integration
points and the results must therefore be processed before writing them into the output file(s).

.. literalinclude:: homPlaneElast.py
   :lines: 110-
