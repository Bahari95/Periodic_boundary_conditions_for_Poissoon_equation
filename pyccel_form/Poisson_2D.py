from simplines import compile_kernel, apply_dirichlet

from simplines import SplineSpace
from simplines import TensorSpace
from simplines import StencilMatrix
from simplines import StencilVector
from simplines import sol_field_2d
from simplines import plot_field_2d

#.. Prologation by knots insertion matrix
from simplines import prolongation_matrix


#from matplotlib.pyplot import plot, show
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
from mpl_toolkits.mplot3d.axes3d import get_test_data
from matplotlib.ticker import LinearLocator, FormatStrFormatter
#%matplotlib inline
import timeit
import time
start = time.time()

from gallery_section_00 import assemble_stiffnessmatrix1D
from gallery_section_00 import assemble_massmatrix1D

assemble_stiffness1D = compile_kernel( assemble_stiffnessmatrix1D, arity=2)
assemble_mass1D      = compile_kernel( assemble_massmatrix1D, arity=2)

#---In Poisson equation
from gallery_section_00 import assemble_vector_ex01 #---1 : In uniform mesh
from gallery_section_00 import assemble_norm_ex01 #---1 : In uniform mesh

assemble_rhs         = compile_kernel(assemble_vector_ex01, arity=1)
assemble_norm_l2  = compile_kernel(assemble_norm_ex01, arity=1)

print('time to import utilities of Poisson equation =', time.time()-start)


from scipy.sparse import kron
from scipy.sparse import csr_matrix
from kronecker.kronecker import vec_2d
from kronecker.fast_diag import Poisson
from scipy.sparse import csc_matrix, linalg as sla
from numpy import zeros, linalg, asarray
import numpy as np

from tabulate import tabulate

#==============================================================================
#.......Poisson ALGORITHM
def poisson_solve(V1, V2 , V, u_p):

       u   = StencilVector(V.vector_space)
       #... We delete the first and the last spline function
       #. as a technic for applying Dirichlet boundary condition
       
       #..Stiffness and Mass matrix in 1D in the first deriction
       K1 = assemble_stiffness1D(V1)
       K1 = K1.tosparse()
       K1 = K1.toarray()
       K1[0,:]  = 0.
       K1[-1,:] = 0.
       #.
       K1[0,0]   = 1.
       K1[-1,-1] = 1.
       K1 = csr_matrix(K1)
       
       M1 = assemble_mass1D(V1)
       M1 = M1.tosparse()
       M1 = M1.toarray()
       M1[0,:]  = 0.
       M1[-1,:] = 0.
       #.
       M1[0,0]   = 1.
       M1[-1,-1] = 1.
       M1 = csr_matrix(M1)

       # Stiffness and Mass matrix in 1D in the second deriction
       K2 = assemble_stiffness1D(V2)
       K2 = K2.tosparse()
       K2 = K2.toarray()
       K2 = csr_matrix(K2)

       M2 = assemble_mass1D(V2)
       M2 = M2.tosparse()
       M2 = M2.toarray()
       M2 = csr_matrix(M2)

       mats_1 = [M1, K1]
       mats_2 = [M2, K2]

       # ...
       poisson = Poisson(mats_1, mats_2)
      
       #---- assemble rhs
       rhs = assemble_rhs( V , fields = [u_p])
       #rhs[0,:]           = 0.
       #rhs[V1.nbasis-1,:] = 0.
       #--Solve a linear system

       b  = rhs.toarray()

       #x  = lu.solve(b)
       #x, inf   = sla.cg(M, b)
       x  = poisson.solve(b)
       x        = x.reshape(V.nbasis)
       x[0, :]  = 0.
       x[-1, :]  = 0.
       u.from_array(V, x)

       Norm    = assemble_norm_l2(V, fields=[u])
       norm    = Norm.toarray()
       l2_norm = norm[0]
       H1_norm = norm[1]
       
       return u, x, l2_norm, H1_norm

degree     = 3
nelements  = 64

#----------------------
# create the spline space for each direction
V1 = SplineSpace(degree=degree, nelements= nelements, nderiv = 2)
V2 = SplineSpace(degree=degree, nelements= nelements, nderiv = 2,periodic=True)

# create the tensor space
Vh = TensorSpace(V1, V2)

# ... Dirichlet boundary condition 
xuh                         = zeros(Vh.nbasis)
xuh[0,:]                    = 0.
u_p                         = StencilVector(Vh.vector_space)
u_p.from_array(Vh, xuh)

print('#---IN-UNIFORM--MESH')
u_ph, xuh, l2_norm, H1_norm = poisson_solve(V1, V2, Vh, u_p)
print('l2_norm = {} H1_norm = {} '.format(l2_norm, H1_norm) ) 
nbpts    = 50

u_poisson, a, b, Y, X = sol_field_2d((nbpts,nbpts),  xuh , Vh.knots, Vh.degree)
figtitle  = 'solution of poisson equation in the new geometry'

fig, axes = plt.subplots( 1, 2, figsize=[12,12], gridspec_kw={'width_ratios': [2.75, 2]} , num=figtitle )
for ax in axes:
   ax.set_aspect('equal')

#axes[0].set_title( 'electric potential Function' )
im = axes[0].contourf( X,Y, u_poisson, cmap= 'jet')
fig.colorbar(im, shrink=0.52, aspect=20)
#axes[1].set_title( 'exact electric potential Function' )
im = axes[1].contourf( X,Y, u_poisson, cmap= 'jet')
fig.colorbar(im, shrink=0.52, aspect=20)
fig.tight_layout()
plt.subplots_adjust(wspace=0.3)
plt.savefig('periodic_Poisson')
plt.show()
