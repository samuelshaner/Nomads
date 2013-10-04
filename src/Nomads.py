
import numpy as np
from math import *
from Cell import *
from Material import *
from Mesh import *
from Plotter import *
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve


# Main solver class. This class sets up the matrix equations
# and solves the finite difference form of the diffusion equation
# using a power method eigenvalue solver 
class Solver(object):

    def __init__(self, mesh):
        
        # initialize variables
        self.mesh = mesh
        self.num_groups = mesh.cells[0].material.num_groups
        self.A = lil_matrix((mesh.num_cells*self.num_groups,mesh.num_cells*self.num_groups))
        self.M = lil_matrix((mesh.num_cells*self.num_groups,mesh.num_cells*self.num_groups))                    
        self.phi = np.ones(mesh.num_cells*self.num_groups)
        self.keff = 1.0
        self.keff_old = 0.0
        

    # Compute the surface diffusion coefficients linking each cell together
    def computeDs(self):
        
        for z in range(self.mesh.cells_z):
            for y in range(self.mesh.cells_y):
                for x in range(self.mesh.cells_x):
                
                    cell = self.mesh.cells[z*self.mesh.cells_x*self.mesh.cells_y + y*self.mesh.cells_x + x]
                
                    for side in range(6):
                    
                        surface = cell.surfaces[side]
                        cellNext = cell.neighborCells[side]
                    
                        for e in range(self.num_groups):
                        
                            d = cell.material.D[e]
                            flux = cell.flux[e]
                        
                            # set the length of the surface parallel to and perpendicular from surface
                            if side == 0 or side == 1:
                                length_perpen = cell.length_x
                            elif side == 2 or side == 3:
                                length_perpen = cell.length_y
                            elif side == 4 or side == 5:
                                length_perpen = cell.length_z
                            
                            # check if cell has a neighbor cell
                            if cellNext is None:
                                # set diffusion coefficients for VACUUM BC
                                if surface.boundary == 'VACUUM':
                                    d_dif = 2 * d / length_perpen / (1 + 4 * d / length_perpen)
                                    d_hat = 2 * d / length_perpen / (1 + 4 * d / length_perpen)
                                # set diffusion coefficients for ZERO FLUX BC
                                elif surface.boundary == 'ZERO_FLUX':
                                    d_dif = 2 * d / length_perpen
                                    d_hat = 2 * d / length_perpen
                                # set diffusion coefficients for REFLECTIVE BC
                                else:
                                    d_hat = 0.0
                                    d_dif = 0.0

                            # if cell has neighbor cell, compute the surface diffusion coefficient                                
                            else:

                                if side == 0 or side == 1:
                                    next_length_perpen = cellNext.length_x
                                elif side == 2 or side == 3:
                                    next_length_perpen = cellNext.length_y
                                elif side == 4 or side == 5:
                                    next_length_perpen = cellNext.length_z

                                d_next = cellNext.material.D[e]
                        
                                d_dif = 2 * d * d_next / (length_perpen * d + next_length_perpen * d_next)
                                d_hat = 2 * d * d_next / (length_perpen * d + next_length_perpen * d_next)
                                                                                            
                            surface.DDif[e] = d_dif
                            surface.DHat[e] = d_hat
                         
    # construct the A (loss) and M (gain) matricies
    def makeAM(self):
        
        cx = self.mesh.cells_x
        cy = self.mesh.cells_y
        cz = self.mesh.cells_z
        ng = self.num_groups

        # zero A and M
        self.A = self.A.tolil()
        self.A = self.A * 0.0
        self.M = self.M.tolil()
        self.M = self.M * 0.0
        
        # loop over mesh cells
        for z in range(cz):
            for y in range(cy): 
                for x in range(cx):
                
                    # loop over energy groups
                    for e in range(ng):
                    
                        cn = z*cy*cx + y*cx + x
                        cell = self.mesh.cells[cn]
                    
                        # absorption term on diagonal
                        self.A[cn*ng + e, cn*ng + e] += cell.material.sigma_a[e] * cell.volume
                    
                        # out scattering term on diagonal
                        for g in range(ng):
                            if e != g:
                                self.A[cn*ng + e, cn*ng + e] += cell.material.sigma_s[e,g] * cell.volume
                    
                        # fission terms on diagonal
                        for g in range(ng):
                            self.M[cn*ng + e, cn*ng + g] = cell.material.chi[e] * cell.material.nu_sigma_f[g] * cell.volume

                        # in scattering terms on off diagonals
                        for g in range(ng):
                            if e != g:
                                self.A[cn*ng + e, cn*ng + g] -= cell.material.sigma_s[g,e] * cell.volume
                            
                        # LEFT SURFACE
                    
                        # transport term on diagonal
                        self.A[cn*ng + e, cn*ng + e] += cell.surfaces[0].DDif[e] * cell.length_y * cell.length_z
                        
                        # transport terms on off diagonals
                        if x != 0:
                            self.A[cn*ng + e, (cn - 1)*ng + e] -= cell.surfaces[0].DDif[e] * cell.length_y * cell.length_z

                        # RIGHT SURFACE
                    
                        # transport term on diagonal
                        self.A[cn*ng + e, cn*ng + e] += cell.surfaces[1].DDif[e] * cell.length_y * cell.length_z
                        
                        # transport terms on off diagonals
                        if x != cx - 1:
                            self.A[cn*ng + e, (cn + 1)*ng + e] -= cell.surfaces[1].DDif[e] * cell.length_y * cell.length_z

                        # BACK SURFACE
                    
                        # transport term on diagonal
                        self.A[cn*ng + e, cn*ng + e] += cell.surfaces[2].DDif[e] * cell.length_x * cell.length_z
                        
                        # transport terms on off diagonals
                        if y != 0:
                            self.A[cn*ng + e, (cn - cx)*ng + e] -= cell.surfaces[2].DDif[e] * cell.length_x * cell.length_z

                        # FRONT SURFACE
                    
                        # transport term on diagonal
                        self.A[cn*ng + e, cn*ng + e] += cell.surfaces[3].DDif[e] * cell.length_x * cell.length_z
                        
                        # transport terms on off diagonals
                        if y != cy - 1:
                            self.A[cn*ng + e, (cn + cx)*ng + e] -= cell.surfaces[3].DDif[e] * cell.length_x * cell.length_z

                        # BOTTOM SURFACE
                    
                        # transport term on diagonal
                        self.A[cn*ng + e, cn*ng + e] += cell.surfaces[4].DDif[e] * cell.length_x * cell.length_y
                        
                        # transport terms on off diagonals
                        if z != cz - 1:
                            self.A[cn*ng + e, (cn + cx*cy)*ng + e] -= cell.surfaces[4].DDif[e] * cell.length_x * cell.length_y

                        # TOP SURFACE
                    
                        # transport term on diagonal
                        self.A[cn*ng + e, cn*ng + e] += cell.surfaces[5].DDif[e] * cell.length_x * cell.length_y
                        
                        # transport terms on off diagonals
                        if z != 0:
                            self.A[cn*ng + e, (cn - cx*cy)*ng + e] -= cell.surfaces[5].DDif[e] * cell.length_x * cell.length_y

        self.A = self.A.tocsr()
        self.M = self.M.tocsr()
        
        
    # solve the eigenvalue problem
    def solve(self, tol):

        # set the BC for each cell
        self.mesh.setCellBoundaries()

        # compute the surface diffusion coefficients
        self.computeDs()
        
        # construct the A (loss) and M (gain) matrices
        self.makeAM()

        # compute the scalar flux
        self.computeFlux(tol)

        print 'DIFFUSION: --- k_eff = ' + str(self.keff)[0:10]


    # solve the generalized eigenvalue problem using the power method
    def computeFlux(self, tol):
                 
        # initialize variables
        max_iter = 3000
        ng = self.mesh.cells[0].material.num_groups
        cx = self.mesh.cells_x
        cy = self.mesh.cells_y
        cz = self.mesh.cells_z
        self.phi = np.ones(cx*cy*cz*ng)
        self.keff_old = self.keff

        # old and new sources
        sold = np.zeros(cx*cy*cz*ng)
        snew = np.zeros(cx*cy*cz*ng)

        # residual
        res = np.zeros(cx*cy*cz*ng)

        # get initial source and find initial keff
        snew = self.M * self.phi
        sumnew = np.sum(snew)
        sold = self.A * self.phi
        sumold = np.sum(sold)
        self.keff = sumnew / sumold
        
        # recompute and initialize the initial source
        sold = self.M * self.phi
        sumold = sum(sold)
        sold = sold * (cx*cy*cz*ng / sumold)
        sumold = cx*cy*cz*ng
        
        # perform power iterations
        for i in range(max_iter):

            # solve Ax=b problem 
            self.phi = spsolve(self.A, sold)
 
            # compute the new source
            snew = self.M * self.phi
            sumnew = np.sum(snew)
            
            # compute and set keff
            self.keff = sumnew / sumold
            sold = sold * self.keff
             
            # compute L2 norm of source error
            snew += 1.e-10
            res = np.divide(sold, snew)
            res = res - 1.0
            l2_norm = np.linalg.norm(res)
            l2_norm = l2_norm / (cx*cy*cz*ng)
         
            # normalize the source and save to old
            snew = snew * (cx*cy*cz*ng) / sumnew             
            sold = snew
            
            print 'iteration: ' + str(i) + '   error: ' + str(abs(self.keff - self.keff_old)) + '  --- keff: ' + str(self.keff)[0:10]
             
            # check for convergence of keff
            if abs(self.keff-self.keff_old) < tol:
                
                # if converged, save the flux
                for i in range(cx*cy*cz):
                    for e in range(ng):
                        self.mesh.cells[i].flux[e] = self.phi[i*ng+e]
                break            

            self.keff_old = self.keff
