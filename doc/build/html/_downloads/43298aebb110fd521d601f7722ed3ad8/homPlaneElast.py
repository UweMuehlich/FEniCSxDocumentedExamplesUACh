# test homogeneous plane elasticity  
# ubuntu Release 18.04.6 LTS (Bionic Beaver) 64-bit
# docker version 24.0.2, build cb74dfc

import numpy as np
import dolfinx.fem.petsc
import ufl
 
from dolfinx  import fem, io, mesh, default_scalar_type
from mpi4py   import MPI
from petsc4py import PETSc

print("===============================================")
print(f"DOLFINx version: {dolfinx.__version__}")
print(f"based on GIT commit: {dolfinx.git_commit_hash}")
print(f"of https://github.com/FEniCS/dolfinx/")
print("===============================================")
#---------------------------------------------------------------------
# DATA SPECIFICATION
#---------------------------------------------------------------------
nr_elemsX1   = 9                   
nr_elemsX2   = 9
strip_height = 1.0 
plate_height = 3.0*strip_height

YoungMod     =  1.e5
PoissonNr    =  0.3
VolumeForce  =  0.0

res_filename = "HomPlaneElastResults.xdmf"
#--------------------------------------------------------------------- 
# MESH
#---------------------------------------------------------------------
my_mesh  = mesh.create_rectangle( comm=MPI.COMM_WORLD,
                                  points=((0.0, 0.0),
                                  (plate_height, plate_height)),
                                  n=(nr_elemsX1, nr_elemsX2),
                                  cell_type= mesh.CellType.quadrilateral)

gdim = my_mesh.geometry.dim

V = fem.functionspace(my_mesh, ("Lagrange", 1, (gdim, )))
#---------------------------------------------------------------------
# SUBDOMAINS DIRICHLET BC'S
#---------------------------------------------------------------------
u_Bottom = np.array([0,     0], dtype=default_scalar_type)
u_Top    = -0.5
u_Left   =  0.0

def Bottom(x):
    return np.isclose(x[1], 0)

def Left(x):
    return np.isclose(x[0], 0)

def Top(x):
    return np.isclose(x[1], plate_height)

fdim = my_mesh.topology.dim - 1

Bottom_facets = mesh.locate_entities_boundary(my_mesh, fdim, Bottom)
Left_facets   = mesh.locate_entities_boundary(my_mesh, fdim, Left)
Top_facets    = mesh.locate_entities_boundary(my_mesh, fdim, Top)

Bottom_dofs = fem.locate_dofs_topological(V, fdim, Bottom_facets)

Left_dofs   = fem.locate_dofs_topological(V.sub(0), fdim, Left_facets)
Top_dofs    = fem.locate_dofs_topological(V.sub(1), fdim, Top_facets)

bcBottom = fem.dirichletbc(  u_Bottom,  Bottom_dofs, V)

bcLeft   = fem.dirichletbc(u_Left, Left_dofs, V.sub(0))
bcTop    = fem.dirichletbc(u_Top , Top_dofs , V.sub(1))

bcs = [bcBottom, bcLeft, bcTop]
#---------------------------------------------------------------------
# CONSTITUTIVE LAW
#---------------------------------------------------------------------
E  = fem.Constant(my_mesh, YoungMod)
nu = fem.Constant(my_mesh, PoissonNr)

model = "plane_stress"

mu_c    = E/2/(1+nu)
lmbda_c = E*nu/(1+nu)/(1-2*nu)

if model == "plane_stress":
    lmbda_c = 2*mu_c*lmbda_c/(lmbda_c+2*mu_c)

def epsilon(u):
    return ufl.sym(ufl.grad(u))

# Stress function
def sigma(u):
    return lmbda_c*ufl.nabla_div(u)*ufl.Identity(2) + 2*mu_c*epsilon(u)

#---------------------------------------------------------------------
# VARIATIONAL PROBLEM
#---------------------------------------------------------------------
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)
f = fem.Constant(my_mesh, default_scalar_type((0, VolumeForce)))
a = ufl.inner(sigma(u), epsilon(v)) * ufl.dx
L = ufl.dot(f, v) * ufl.dx 

petsc_opts={"ksp_type": "preonly", "pc_type": "lu"}

problem = fem.petsc.LinearProblem(a, L, bcs, petsc_options=petsc_opts)
uh      = problem.solve()
#---------------------------------------------------------------------
# RESULTS
#---------------------------------------------------------------------
V0      = fem.functionspace(my_mesh, ("Lagrange", 1, (gdim, gdim)))
sig_exp = fem.Expression(sigma(uh), V0.element.interpolation_points())
sig     = fem.Function(V0, name="Stress")
sig.interpolate(sig_exp)

out_ASCII = False
if out_ASCII:
   encoding = dolfinx.io.XDMFFile.Encoding.ASCII
else: 
   encoding = dolfinx.io.XDMFFile.Encoding.HDF5


xdmffile = dolfinx.io.XDMFFile(MPI.COMM_WORLD, res_filename, "w", encoding)

xdmffile.write_mesh(my_mesh)
uh.name = "Displacement"
xdmffile.write_function(uh)
#sig.name = "Stress"
#xdmffile.write_function(sig)
xdmffile.close()



