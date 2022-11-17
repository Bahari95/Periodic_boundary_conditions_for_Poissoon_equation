from evtk.hl import gridToVTK 
import numpy as np 
import random as rnd 

# Coordinates
x = np.load("mesh_x.npy")
y = np.load("mesh_y.npy")
z = np.load("mesh_z.npy")
# We add Jacobian function to make the grid more interesting
Ex_sol       = np.load("Ex_sol.npy")
appr_sol     = np.load("appr_sol.npy")
appr_soldx   = np.load("appr_soldx.npy")
appr_soldy   = np.load("appr_soldy.npy")
appr_soldz   = np.load("appr_soldz.npy")
# Dimensions 
nx, ny, nz = x.shape[0], y.shape[0], z.shape[0]
ncells = nx * ny * nz
npoints = (nx + 1) * (ny + 1) * (nz + 1)

# Variables <br>
gridToVTK("./domain", x, y, z, pointData = {"Ex_sol" : Ex_sol, "appr_sol" : appr_sol,"dx_appr_sol" : appr_soldx,"dy_appr_sol" : appr_soldy,"dz_appr_sol" : appr_soldz,})
