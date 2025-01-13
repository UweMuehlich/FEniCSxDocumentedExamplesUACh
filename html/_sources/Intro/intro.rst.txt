Introduction
=============
The Pros and Cons of Using FEniCSx for Research in Engineering
--------------------------------------------------------------

FEniCS has long been a well-regarded tool in the research community for solving
partial differential equations (PDEs). With the development of FEniCSx, there is
a strong expectation that it will be embraced in the same way, if not more so.

The FEniCSx open-source project reflects the work of highly skilled and ambitious
developers, and its integration into contemporary research is a testament to
their efforts

Advantages of using FEniCSx
...........................

One of the most appealing aspects of FEniCSx is its high-level programming
interface, particularly the Unified Form Language (UFL), which simplifies the
process of formulating and solving PDEs. For standard problems, the ability to
create solutions with just a few lines of intuitive code is impressive.

This efficiency allows researchers to focus more on the problem at hand rather
than getting bogged down in implementation details. The ease with which you can
solve relatively simple problems with FEniCSx is one of its key strengths.

A key advantage of FEniCSx is its open-source nature, offering a free
alternative to costly commercial FEM software. While commercial tools are seen
as mature but expensive, FEniCSx enables academic researchers to collaborate
and innovate without financial constraints.

Challenges of using FEniCSx
...........................

While FEniCSx shines in handling standard problems, its use becomes more
challenging when addressing non-standard or highly specialized research
questions. To effectively apply FEniCSx to these types of problems, a deep
understanding is required—not only of the theoretical foundations of FEM but
also of the software’s inner workings.
This dual requirement can present a barrier for students and researchers
who are not computational experts.

Planned Topics and Examples
---------------------------

.. note::
   This page is a work in progress and may be updated occasionally.

Linear Elastostatics

* :ref:`Basic plane problem <Target_MinimalExamplePES>`:   ✔
* :ref:`Plane problem with two isotropic materials <Target_PES2Mat>`:   ✔
* :ref:`Plane problem with two orthotropic materials <Target_PESOrtho2Mat>`:   ✔
* Basic 3-d problem
* 3-d problem with different materials
* Homogenization method with fluctuations as primary unknowns (method I)
* Homogenization method with total displacements as primary unknowns (method II)
  using the multiple-point-constraint feature of Dolfinx

Linear Elastodynamics

* :ref:`Modal analysis, minimal example <Target_ModalGA>`:   ✔
* Modal analysis in-plane
* Modal analysis in-plane with hole
* Dispersion curves 1-d periodic cells
* Dispersion curves 2-d periodic cells
* Dispersion curves 3-d periodic cells
* Modal analysis CLT plate
