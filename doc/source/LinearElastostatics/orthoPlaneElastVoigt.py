# Orthotropic plane elasticity  Voigt Notation
# ubuntu Release 18.04.6 LTS (Bionic Beaver) 64-bit
# docker version 24.0.2, build cb74dfc

import numpy as np
import dolfinx.fem.petsc
import ufl
 
from dolfinx  import fem, io, mesh, default_scalar_type
from mpi4py   import MPI
from petsc4py import PETSc

folder        = "../Results/"
problem_label = "OrthoPlaneElast2Mats_"

print("===============================================")
print(f"DOLFINx version: {dolfinx.__version__}")
print(f"based on GIT commit: {dolfinx.git_commit_hash}")
print(f"of https://github.com/FEniCS/dolfinx/")
print(f"Problem label: {problem_label}")
print("===============================================")
#---------------------------------------------------------------------
# DATA SPECIFICATION
#---------------------------------------------------------------------
nr_elemsX1    = 200                
nr_elemsX2    = 50
strip_height  = 1.0 
plate_height  = 5.0*strip_height
plate_width   = 10.0*plate_height

E2, E1, nu12, G12 = 7000., 300., 0.38, 260.

VolumeForce  =  0.0

res_filename  = folder + problem_label + "Results.xdmf"
#--------------------------------------------------------------------- 
# MESH
#---------------------------------------------------------------------
my_mesh  = mesh.create_rectangle( comm=MPI.COMM_WORLD,
                                  points=((0.0, 0.0),
                                  (plate_width, plate_height)),
                                  n=(nr_elemsX1, nr_elemsX2),
                                  cell_type= mesh.CellType.quadrilateral)

gdim = my_mesh.geometry.dim

V = fem.functionspace(my_mesh, ("Lagrange", 1, (gdim, )))
#---------------------------------------------------------------------
# DIRICHLET BC'S
#---------------------------------------------------------------------
u_Bottom = np.array([0,     0], dtype=default_scalar_type)
u_Top    = -0.05*plate_height
u_Left   =  np.array([0,     0], dtype=default_scalar_type)

def Left(x):
    return np.isclose(x[0], 0)

def Right(x):
    return np.isclose(x[0], plate_width)

fdim = my_mesh.topology.dim - 1

Left_facets   = mesh.locate_entities_boundary(my_mesh, fdim, Left)
Right_facets  = mesh.locate_entities_boundary(my_mesh, fdim, Right)


Left_dofs   = fem.locate_dofs_topological(V, fdim, Left_facets)
Right_dofs  = fem.locate_dofs_topological(V.sub(1), fdim, Right_facets)

bcLeft   = fem.dirichletbc(u_Left, Left_dofs, V)
bcRight  = fem.dirichletbc(u_Top , Right_dofs, V.sub(1))

bcs = [bcLeft, bcRight]
#---------------------------------------------------------------------
# DOMAINS FOR DIFFERENT MATERIALS
#---------------------------------------------------------------------
def EvenLayer(x):
    in_Even   = [ int(y // strip_height ) % 2 == 0     for y in x[1]]
    at_Border = [ np.isclose(y, int(y / strip_height)) for y in x[1]]
    return np.array(in_Even) | np.array(at_Border) 

cells_Even = mesh.locate_entities(my_mesh, my_mesh.topology.dim, EvenLayer)

# Get the total number of cells in the mesh
num_cells = my_mesh.topology.index_map(my_mesh.topology.dim).size_local

# Get all cell indices (0 to num_cells-1)
all_cells = np.arange(num_cells)

Q = fem.functionspace(my_mesh, ("DG", 0))

indicator = fem.Function(Q)

indicator.x.array[all_cells ] = np.full_like(all_cells , 0.0, dtype=PETSc.ScalarType)
indicator.x.array[cells_Even] = np.full_like(cells_Even, 1.0, dtype=PETSc.ScalarType)
#---------------------------------------------------------------------
# CONSTITUTIVE LAW
#---------------------------------------------------------------------
S  = np.array([[1./E1,-nu12/E1,0.],[-nu12/E1,1./E2,0.],[0.,0.,1./G12]])
ST = np.array([[1./E2,-nu12/E2,0.],[-nu12/E2,1./E1,0.],[0.,0.,1./G12]])

C  = np.linalg.inv(S)
CT = np.linalg.inv(ST)

AA = PETSc.ScalarType( ( (C[0,0], C[0,1], C[0,2]), 
                         (C[1,0], C[1,1], C[1,2]), 
                         (C[2,0], C[2,1], C[2,2])  ) )

BB = PETSc.ScalarType( ( (CT[0,0], CT[0,1], CT[0,2]), 
                         (CT[1,0], CT[1,1], CT[1,2]), 
                         (CT[2,0], CT[2,1], CT[2,2])  ) )

CC0 = fem.Constant( my_mesh,AA ) 
CC1 = fem.Constant( my_mesh,BB ) 

def strain(u, repr="vectorial"):
    eps_t = ufl.sym(ufl.grad(u))
    if repr == "vectorial":
        return ufl.as_vector([eps_t[0, 0], eps_t[1, 1], 2 * eps_t[0, 1]])
    elif repr == "tensorial":
        return eps_t

def stress(u, repr="vectorial"):
    sigv = (1.0-indicator)*ufl.dot(CC0,strain(u)) + indicator*ufl.dot(CC1,strain(u))
    if repr == "vectorial":
        return sigv
    elif repr == "tensorial":
        return ufl.as_matrix([[sigv[0], sigv[2]], [sigv[2], sigv[1]]])
#---------------------------------------------------------------------
# VARIATIONAL PROBLEM
#---------------------------------------------------------------------
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)
f = fem.Constant(my_mesh, default_scalar_type((0, VolumeForce)))
a = ufl.inner(stress(u), strain(v)) * ufl.dx
L = ufl.dot(f, v) * ufl.dx 

petsc_opts={"ksp_type": "preonly", "pc_type": "lu"}

problem = fem.petsc.LinearProblem(a, L, bcs, petsc_options=petsc_opts)
uh      = problem.solve()
#---------------------------------------------------------------------
# OUTPUT RESULTS
#---------------------------------------------------------------------
V0      = fem.functionspace(my_mesh, ("Lagrange", 1, (gdim, gdim)))
sig_exp = fem.Expression(stress(uh,"tensorial"), V0.element.interpolation_points())
sig     = fem.Function(V0, name="Stress")
sig.interpolate(sig_exp)

eps_exp = fem.Expression(strain(uh,"tensorial"), V0.element.interpolation_points())
eps     = fem.Function(V0, name="Strain")
eps.interpolate(eps_exp)


vtkD = io.VTKFile(my_mesh.comm, folder+ problem_label + "Displacements.pvd", "w")
vtkD.write_function(uh)
vtkD.close()

vtkS = io.VTKFile(my_mesh.comm, folder + problem_label + "Stress.pvd", "w")
vtkS.write_function(sig)
vtkS.close()

vtkE = io.VTKFile(my_mesh.comm, folder + problem_label + "Strain.pvd", "w")
vtkE.write_function(eps)
vtkE.close()

vtkM = io.VTKFile(my_mesh.comm, folder + problem_label + "Materials.pvd", "w")
vtkM.write_function(indicator)
vtkM.close()


