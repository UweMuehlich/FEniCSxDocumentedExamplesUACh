Example using SLEPc  for a complex eigenvalue problem
------------------------------------------------------------

We consider a beam-like structure clamped at one end, assuming a linear elastic
material. For illustration purpose, we try to solve the problem

.. math::
   \begin{eqnarray}
     \left[K + \alpha \mathrm{I} B \right] \mathbf{x} = \omega^2 M \mathbf{x}
   \end{eqnarray}

with  :math:`B=K`, where :math:`\mathrm{I}` is the imaginary unit.

.. note::
   You can download the complete python script
   :download:`TestModalComplex.py <codeGA/TestModalComplex.py>`

**Imports and Environment Check**

.. literalinclude:: codeGA/TestModalComplex.py
   :lines: 8-22

* Imports essential libraries for mesh creation, variational formulation, and eigenvalue solving.
* Ensures complex arithmetic (numpy.complex128) is en\mathrm{I}abled, as required for this problem.

**Function to Create Complex Stiffness Matrix**

.. literalinclude:: codeGA/TestModalComplex.py
   :lines: 23-52

* Avoids MatCopy Issue: Instead of copying K and modifying it, create_complex_stiffness directly
  constructs the matrix using setValuesCSR, ensuring compatibility with PETSc's complex arithmetic.
* Checks Matrix Dimensions: Ensures K and B are compatible.


**Domain and Mesh Definition**

.. literalinclude:: codeGA/TestModalComplex.py
   :lines: 53-59

Defines a 3D hexahedral mesh for a beam-like domain.

**Constitutive Law**

.. literalinclude:: codeGA/TestModalComplex.py
   :lines: 60-72

Defines the material properties and constitutive relations for linear elasticity,
assuring compatibility with complex arithmetic using PETSc.ScalarType.

**Function Space and Boundary Conditions**

.. literalinclude:: codeGA/TestModalComplex.py
   :lines: 73-86

Defines the function space, trial and test functions, and clamped boundary conditions.

**Assemble Stiffness and Mass Matricess and Complex Matrix A**

.. literalinclude:: codeGA/TestModalComplex.py
   :lines: 87-104

Assembles the stiffness and mass matrices using variational forms
and creates a complex-valued stiffness matrix
using the previously defined function.

**Solve the Eigenvalue Problem and Writes the Results**

.. literalinclude:: codeGA/TestModalComplex.py
   :lines: 105-


**Remarks**

* Direct matrix operations (e.g., MatCopy) may fail if the matrix is not in a
  fully compatible state. This was resolved by using setValuesCSR to explicitly
  construct the complex matrix.
* Enabling complex128 mode ensures compatibility for problems involving
  imaginary units.

.. note::
   In case of an error  message try first **source dolfinx-complex-mode** in your
   Docker environment.
