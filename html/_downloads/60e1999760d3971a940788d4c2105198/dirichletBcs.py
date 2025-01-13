# DOLFINX tutorial, topic Dirichlet bcs
"""
A box domain defined by [0,L0]x[0,L1]x[0,L2] is considered. 

Linear elasticity is supposed in the interior of the box.

The following Dirichlet boundary conditions will be prescribed:
 I)   The facet {0}x[0,L1]x[0,L2] (facets on the x0 = 0 plane) is completely constrained,
      u  = 0 for X in  {0}x[0,L1]x[0,L2] where u =[u0, u1, u2]
 II)  u3 = 0 for X in {L0}x[0,L1]x{0}    (edges on the line (x0=L0, x1, x2=0)
 III) [u2=v, u3 = v] for X = {L0}x{L2}x{L2} (corner point with x = (L0,L1,L2)
Note: X coordinate, x Cartesian product
"""
import numpy as np
import ufl
import dolfinx.fem.petsc

from dolfinx import fem, mesh, default_scalar_type
from mpi4py  import MPI

print("--------------------------------------------------------------")
print("--------    ILLUSTRATING DIRICHLET BC'S   --------------------")
print("--------------------------------------------------------------")

L0, L1, L2 = 10., 5., 1.
n0, n1, n2 = 10,  5,  4
my_mesh  = mesh.create_box( comm=MPI.COMM_WORLD,                    # creating mesh
                            points=((0.0, 0.0, 0.0),
                                    (L0, L1, L2)),
                            n=(n0, n1, n2),
                            cell_type= mesh.CellType.hexahedron)

gdim = my_mesh.geometry.dim
tdim = my_mesh.geometry.dim
print(f"gdim {gdim}, t_dim {tdim}")

V = fem.functionspace(my_mesh, ("Lagrange", 1, (gdim, )))           # creating  
                                                                    # function space
print(f"V.element {type(V.element)}")
print(f"V.dofmap        {type(V.dofmap)}")

#--------------------------------------------------------------------
# BOUNDARY CONDITIONS
#
# I) block all displacements on the x[0]=0 plane
#--------------------------------------------------------------------
def True_if_x0_is_zero(x):
    """ locate if on the x0 = 0.0 plane """
    #print(x.T) # turns out to be just to be a node-coordinate list
    #print(np.isclose(x[0], 0.0))
    return  np.isclose(x[0], 0.0)

dof_g_list = fem.locate_dofs_geometrical(V, True_if_x0_is_zero)  # finds all points but
                                                                  # cannot work with
                                                                  # subspace V.sub(i)

boundary_facets = mesh.locate_entities_boundary(my_mesh, gdim-1, True_if_x0_is_zero)
boundary_edges  = mesh.locate_entities_boundary(my_mesh, gdim-2, True_if_x0_is_zero)
boundary_nodes  = mesh.locate_entities_boundary(my_mesh, gdim-3, True_if_x0_is_zero)

print(f"Dirichlet I")  
print(f"boundary facets on X[0] = 0 plane:  {boundary_facets}")  
print(f"boundary edges  on X[0] = 0 plane:  {boundary_edges}")
print(f"boundary nodes  on X[0] = 0 plane:  {boundary_nodes}")

facets_dofs = fem.locate_dofs_topological(V, gdim-1, boundary_facets)

my_mesh.topology.create_connectivity(gdim-2, gdim) 

edge_dofs   = fem.locate_dofs_topological(V, gdim-2, boundary_edges)

my_mesh.topology.create_connectivity(gdim-3, gdim) 

node_dofs   = fem.locate_dofs_topological(V, gdim-3, boundary_nodes)

print("\nAll four methods give the same dofs")
print(f"geom_dofs   on X[0] = 0 plane:  {dof_g_list}")  
print(f"facets_dofs on X[0] = 0 plane:  {facets_dofs}")  
print(f"edge_  dofs on X[0] = 0 plane:  {edge_dofs}")
print(f"node  _dofs on X[0] = 0 plane:  {node_dofs}")

u_bcI    = np.array([0, 0, 0], dtype=default_scalar_type)

# The following four methods lead the same result
#bcI      = fem.dirichletbc(u_bcI,dof_g_list, V)
#bcI      = fem.dirichletbc(u_bcI,bcI_dofs, V)
#bcI      = fem.dirichletbc(u_bcI,dof_g_list, V)
bcI       = fem.dirichletbc(u_bcI,dof_g_list, V)
#--------------------------------------------------------------------
# II) prescribe u[2]= 0.0  at the edge with coordinates (L0, x[1], 0)
#--------------------------------------------------------------------
def True_if_on_edge(x):
    """ locates if on the edge with (L0, x[1], 0) """
    return np.logical_and(np.isclose(x[0], L0), np.isclose(x[2], 0.0))

edge_nodes  = mesh.locate_entities_boundary(my_mesh, gdim-2, True_if_on_edge)

my_mesh.topology.create_connectivity(gdim-2, gdim) 

edge_dofs  = fem.locate_dofs_topological(V, gdim-2, edge_nodes)
edge_dofs2 = fem.locate_dofs_topological(V.sub(2), gdim-2, edge_nodes)

print(f"\nDirichlet II")  
print(f"edge_dofs  {edge_dofs}")
print(f"edge_dofs2 {edge_dofs2}")

u_bc_edge   = 0.0
bc_edge     = fem.dirichletbc(u_bc_edge,edge_dofs2, V.sub(2))
#--------------------------------------------------------------------
# III) prescribe u[0], u[1] and u[2] at corner node (L0,L1,L2)
#--------------------------------------------------------------------
def True_if_corner(x):
    """ locates if corner point (L0, L1, L2) """
    check = np.logical_and(
              np.logical_and(np.isclose(x[0], L0), np.isclose(x[1], L1)), 
              np.isclose(x[2], L2))
    #print(f"check {check}\n")
    return(check)

corner_node  = mesh.locate_entities_boundary(my_mesh, gdim-3, True_if_corner)

dg_corner =  fem.locate_dofs_geometrical(V, True_if_corner) 

corner_dof0 = fem.locate_dofs_topological(V.sub(0), gdim-3, dg_corner)
corner_dof1 = fem.locate_dofs_topological(V.sub(1), gdim-3, dg_corner)
corner_dof2 = fem.locate_dofs_topological(V.sub(2), gdim-3, dg_corner)

print(f"\nDirichlet III")  
print(f"\ncorner_node {corner_node}")
print(f"dg_corner   {dg_corner}")
print(f"corner_dof0 {corner_dof0}")
print(f"corner_dof1 {corner_dof1}")
print(f"corner_dof2 {corner_dof2}")

u0_corner = 1.0*L2
u1_corner = 1.0*L2
u2_corner = 1.0*L2

bc_corner_0      = fem.dirichletbc(u0_corner,corner_dof0, V.sub(0))
bc_corner_1      = fem.dirichletbc(u1_corner,corner_dof1, V.sub(1))
bc_corner_2      = fem.dirichletbc(u2_corner,corner_dof2, V.sub(2))
#--------------------------------------------------------------------
#    CONSTITUTIVE LAW
#--------------------------------------------------------------------
mu      = 1000.
lamda   = 1000.
rho     = 1.0
gravity = 9.0 

def epsilon(u):
    return ufl.sym(ufl.grad(u))  

def sigma(u):
    return lamda * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2 * mu * epsilon(u)
#--------------------------------------------------------------------
#   VARIATIONAL PROBLEM
#--------------------------------------------------------------------
u  = ufl.TrialFunction(V)
v  = ufl.TestFunction(V)
ff = fem.Constant(my_mesh, default_scalar_type((0, 0, -rho*gravity)))
a  = ufl.inner(sigma(u), epsilon(v)) * ufl.dx
L  = ufl.dot(ff, v) * ufl.dx 

petsc_opts={"ksp_type": "preonly", "pc_type": "lu"}

bcs = [bcI , bc_edge, bc_corner_0, bc_corner_1,  bc_corner_2]
problem = fem.petsc.LinearProblem(a, L, bcs=bcs, petsc_options=petsc_opts)
uh = problem.solve()
uh.name = "Displacement"
#--------------------------------------------------------------------
#   OUTPUT FOR PARAVIEW
#--------------------------------------------------------------------
xdmffile = dolfinx.io.XDMFFile(MPI.COMM_WORLD, "Results.xdmf", "w")

xdmffile.write_mesh(my_mesh)
xdmffile.write_function(uh)
xdmffile.close()

