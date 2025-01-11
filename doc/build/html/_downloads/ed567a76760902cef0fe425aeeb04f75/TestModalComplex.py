"""
 USAGE OF SciPy EIGENSOLVER FOR A PROBLEM
(K+I*B)x=omega*M WITH IMAGINARY UNIT I.
"""
#--------------------------------------------------------------------
# IMPORTS
#--------------------------------------------------------------------
import numpy as np
import ufl
from dolfinx import fem, mesh, default_scalar_type
from dolfinx.fem.petsc import assemble_matrix
from mpi4py import MPI
from petsc4py import PETSc
from slepc4py import SLEPc

print("================START===================")
print(PETSc.ScalarType)
print("If not <class 'numpy.complex128'> run:")
print("source dolfinx-complex-mode") 
print("in your Docker environment")
print("========================================")
#--------------------------------------------------------------------
# FUNCTION TO CREATE COMPLEX MATRIX A = K+j*B (j imag. unit)
#--------------------------------------------------------------------    
def create_complex_stiffness(K,B):
    """Create a complex stiffness matrix A = K + j*B."""
    # Ensure K and B are fully assembled
    if not K.assembled:
        K.assemble()
    if not B.assembled:
        B.assemble()
    # Check if  K and B are of same size
    if not K.getSize()== B.getSize(): 
       print("Error in create_complex_stiffness: K and B must be of same size")

    # Create a new PETSc matrix for the complex stiffness
    A = PETSc.Mat().createAIJ(size=K.getSize(), comm=MPI.COMM_WORLD)
    A.setUp()

    # Get the CSR representation of K
    rows, cols = K.getSize()
    ai, aj, av = K.getValuesCSR()
    bi, bj, bv = K.getValuesCSR()

    # Create a complex array for the new values
    complex_values = av + 1j * bv

    # Assemble A with the complex values
    A.setValuesCSR(ai, aj, complex_values)
    A.assemble()
    
    return A
#--------------------------------------------------------------------
# DOMAIN AND MESH
#--------------------------------------------------------------------    
L = [5.0, 0.6, 0.4]
N = [25, 3, 2]
beam_mesh = mesh.create_box(MPI.COMM_WORLD, [np.zeros(3), np.array(L)], 
                            N, mesh.CellType.hexahedron)
#--------------------------------------------------------------------
# CONSTITUTIVE LAW (LINEAR ELASTICITY)
#--------------------------------------------------------------------
E, nu, rho = 2e11, 0.3, 7850
# Define complex-valued constants for linear elasticity
mu    = fem.Constant(beam_mesh, PETSc.ScalarType(E / (2 * (1 + nu))))
lamda = fem.Constant(beam_mesh, PETSc.ScalarType(E * nu / ((1 + nu) * (1 - 2 * nu))))

def epsilon(u):
    return ufl.sym(ufl.grad(u))

def sigma(u):
    return lamda * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2 * mu * epsilon(u)
#--------------------------------------------------------------------
# FUNCTION SPACE, TRIAL AND TEST FUNCTION
#--------------------------------------------------------------------
V = fem.functionspace(beam_mesh, ("CG", 1,(beam_mesh.geometry.dim, )))
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)
#--------------------------------------------------------------------
# DIRICHLET BOUNDARY CONDITIONS
#--------------------------------------------------------------------
def clamped_boundary(x): return np.isclose(x[0], 0)

u_zero = np.array((0,) * beam_mesh.geometry.dim, dtype=default_scalar_type)

bc = fem.dirichletbc(u_zero, fem.locate_dofs_geometrical(V, clamped_boundary), V)
#--------------------------------------------------------------------
# VARIATIONAL FORMS, STIFFNESS AND MASS MATRIX
#--------------------------------------------------------------------
k_form = fem.form(ufl.inner(sigma(u), epsilon(v)) * ufl.dx)
m_form = fem.form(rho * ufl.inner(u, v) * ufl.dx)

# Assemble matrices
K = assemble_matrix(k_form, bcs=[bc], diagonal=62831)
K.assemble()

M = assemble_matrix(m_form, bcs=[bc], diagonal=1./62831)
M.assemble()
#--------------------------------------------------------------------
# CREATE COMPLEX STIFFNESS MATRIX
#--------------------------------------------------------------------

A = create_complex_stiffness(K,K)

#--------------------------------------------------------------------
# SLEPc EIGENSOLVER 
#--------------------------------------------------------------------
#
# Set up SLEPc eigenvalue solver
eigensolver = SLEPc.EPS().create(MPI.COMM_WORLD)
eigensolver.setOperators(A, M)
eigensolver.setProblemType(SLEPc.EPS.ProblemType.GHEP)
eigensolver.setDimensions(16)  # Number of eigenpairs to compute
eigensolver.setFromOptions()

# Solve the eigenvalue problem
eigensolver.solve()
#--------------------------------------------------------------------
# OUTPUT RESULTS
#--------------------------------------------------------------------
n_converged = eigensolver.getConverged()
if n_converged > 0:
    vr, vi = A.getVecs()
    for i in range(min(n_converged, 16)):
        eigval = eigensolver.getEigenpair(i, vr, vi)
        omega = np.sqrt(eigval.real)
        print(f"Eigenvalue {i}: λ = {eigval:.6e}, ω = {omega:.6f} rad/s")

