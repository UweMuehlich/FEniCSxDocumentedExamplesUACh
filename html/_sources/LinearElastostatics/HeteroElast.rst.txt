.. _Target_PES2Mat:

Plane elastostatics with two different isotropic materials
----------------------------------------------------------

We depart from the :ref:`minimal example for plane elastostatics <Target_MinimalExamplePES>`.
The domain is divided into three horizontal stripes, where the elastic properties
of the central stripe differ from the remaining part.

.. note::
   You can download the complete python script
   :download:`heteroPlaneElast.py <heteroPlaneElast.py>`

Apart from changing the name of the output file, only the section **Constitutive Law**
is affected. Therefore, we discuss in detail only this section.

.. literalinclude:: heteroPlaneElast.py
   :lines: 76-126

The visualization produced by paraview is shown below.

.. figure:: heteroElast.png
   :width: 600

   Deformed configuration, where the colors indicate the displacement in
   vertical direction.
