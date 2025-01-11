.. _Target_PESOrtho2Mat:

Plane Elastostatics with Two Different Orthotropic Materials
------------------------------------------------------------

We adjust the :ref:`minimal example for plane elastostatics <Target_MinimalExamplePES>`.
The domain is modified to a beam-like geometry and divided into horizontal layers (stripes). An orthotropic material is considered, with the principal directions alternating between the layers.

The beam is clamped at the left hand side and for all nodes at the
right hand side a vertical displacement is prescribed.

The visualization produced by paraview is shown below.

.. figure:: orthoPlaneElast.png
   :width: 800

   The following data are shown: deformed configuration together with the displacement magnitude (upper left), the indicator function
   to distinguish between different layers (upper right), stress
   :math:`\sigma_{11}`  (lower left), and strain   :math:`e_{11}`
   (lower right).

.. note::
   The example for FEniCS published by `Jeremy Bleyer <https://comet-fenics.readthedocs.io/en/latest/demo/elasticity/orthotropic_elasticity.py.html>`_
   served as an inspiration for this problem.

   You can download the complete python script
   :download:`orthoPlaneElastVoigt.py <orthoPlaneElastVoigt.py>`




Only the most decisive parts for defining layers of different orthotropic materials
are discussed in detail in the following.

**Defining material layers**

.. literalinclude:: orthoPlaneElastVoigt.py
   :lines: 74-95

We are working with two different materials using an indicator function defined  on a piecewise constant functions space.
If a cell belongs to an even layer (0, 2, ...), the indicator function takes the value one, otherwise zero. In order to identify correctly the cells belonging to a layer, it is important to incoporate positions at the layer boundaries in the corresponding function (here, EvenLayer).


**Constitutive relations**

.. literalinclude:: orthoPlaneElastVoigt.py
   :lines: 96-128


Here, we are using Voigt-notation, defining the compliance of the material with respect to the global axes, :math:`x_1` and :math:`x_2`, using numpy arrays. The compliance matrix for the odd layers is given by

.. math::
      \begin{eqnarray}
      \begin{bmatrix} e_{11} \\ e_{22} \\ 2 e_{12} \end{bmatrix} &=&
      \begin{bmatrix} 1/E_1 & -\frac{\nu_{12}}{E1} & 0 \\
                     -\frac{\nu_{12}}{E1} & 1/E_1   & 0 \\
                     0 & 0 & \frac{1}{G_{12}}
      \end{bmatrix}
      =  \begin{bmatrix} \sigma_{11} \\ \sigma_{22} \\  \sigma_{12} \end{bmatrix}
      \end{eqnarray}

whereas for the even layers, we just switch :math:`E_1` with :math:`E_2`. It means, that material in the even layers are rotated by 90 degree with respect to the :math:`x_3` axis.  After inverting the compliance matrices to compute the corresponding stiffness matrices, we construct the PETSc matrices. By leveraging the indicator function to compute stresses, the appropriate material properties are assigned to each layer.

**Variational problem and result output**

.. literalinclude:: orthoPlaneElastVoigt.py
   :lines: 129-
