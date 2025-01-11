Complex Eigenvalue Problems
---------------------------

In the context of computing dispersion curves for waveguides, eigenvalue problems
of the type

.. math::
   \begin{eqnarray}
     \left[K + \alpha \mathrm{I} B \right] \mathbf{x} = \omega^2 M \mathbf{x}
   \end{eqnarray}

have to be solved, where  :math:`K, B, M` are real Hermitian matrices, :math:`\alpha` is some real number,
and :math:`\mathrm{I}` is the imaginary unit.

There are at least two ways to proceed. In line with the general philosophy of FEniCSx, leveraging the capabilities of SLEPc for handling complex numbers would be the most consistent choice, particularly for large-scale problems. However, this approach can be somewhat challenging initially. An alternative, more straightforward option suitable for small to medium-sized problems, is to use SciPy. In this discussion, we explore both methods.

Before we start, explaining involved acronyms might be helpful:

* PETSc (Portable, Extensible Toolkit for Scientific Computation)
* SLEPc (Scalable Library for Eigenvalue Problem Computations)
* SciPy CSR Matrices (Compressed Sparse Row Format)

.. toctree::
   :maxdepth: 2

   exampleSciPyReal
   exampleSLEPcComplex

In your docker environment for dolfinx/dolfinx:stable, you have to run the command source dolfinx-complex-mode.
