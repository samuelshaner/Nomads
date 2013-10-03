
import numpy as np
import matplotlib.pyplot as plt
from math import *
from Cell import *
from Material import *
from Mesh import *
from Plotter import *
import os
import sys
import getopt
import time
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve


class Solver(object):

    def __init__(self, mesh):
        
        self.mesh = mesh
        self.ng = mesh.cells[0].material.num_groups
        self.A = lil_matrix((mesh.num_cells*self.ng,mesh.num_cells*self.ng))
        self.M = lil_matrix((mesh.num_cells*self.ng,mesh.num_cells*self.ng))
                    
        self.phi = np.ones(mesh.num_cells*self.ng)
        self.keff = 1.0
        self.keff_old = 0.0
        self.relax = 1.0
        
    def computeDs(self):
        
        cx = self.mesh.cells_x
        cy = self.mesh.cells_y
        cz = self.mesh.cells_z

        flux = 0.0 

        for z in range(cz):
            for y in range(cy):
                for x in range(cx):
                
                    cell = self.mesh.cells[z*cx*cy + y*cx + x]
                
                    for side in range(6):
                    
                        surface = cell.surfaces[side]
                        cellNext = cell.neighborCells[side]
                    
                        for e in range(self.ng):
                        
                            d = cell.material.D[e]
                            flux = cell.flux[e]
                        
                            # set the length of the surface parallel to and perpendicular from surface
                            if side == 0 or side == 1:
                                length_perpen = cell.width
                            elif side == 2 or side == 3:
                                length_perpen = cell.length
                            elif side == 4 or side == 5:
                                length_perpen = cell.height
                            
                            if cellNext is None:
                                if surface.boundary == 'VACUUM':
                                    d_dif = 2 * d / length_perpen / (1 + 4 * d / length_perpen)
                                    d_hat = 2 * d / length_perpen / (1 + 4 * d / length_perpen)
                                elif surface.boundary == 'ZERO_FLUX':
                                    d_dif = 2 * d / length_perpen
                                    d_hat = 2 * d / length_perpen
                                else:
                                    d_hat = 0.0
                                    d_dif = 0.0
                                
                            else:

                                if side == 0 or side == 1:
                                    next_length_perpen = cellNext.width
                                elif side == 2 or side == 3:
                                    next_length_perpen = cellNext.length
                                elif side == 4 or side == 5:
                                    next_length_perpen = cellNext.height

                                d_next = cellNext.material.D[e]
                        
                                d_dif = 2 * d * d_next / (length_perpen * d + next_length_perpen * d_next)
                                d_hat = 2 * d * d_next / (length_perpen * d + next_length_perpen * d_next)
                                                                                            
                            surface.DDif[e] = d_dif
                            surface.DHat[e] = d_hat
                         
                         
    def makeAM(self):
        
        cx = self.mesh.cells_x
        cy = self.mesh.cells_y
        cz = self.mesh.cells_z
        ng = self.ng

        self.A = self.A.tolil()
        self.A = self.A * 0.0
        self.M = self.M.tolil()
        self.M = self.M * 0.0
        
        for z in range(cz):
            for y in range(cy): 
                for x in range(cx):
                
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
                            
                        # RIGHT SURFACE
                    
                        # transport term on diagonal
                        self.A[cn*ng + e, cn*ng + e] += cell.surfaces[1].DDif[e] * cell.length * cell.height
                        
                        # transport terms on off diagonals
                        if x != cx - 1:
                            self.A[cn*ng + e, (cn + 1)*ng + e] -= cell.surfaces[1].DDif[e] * cell.length * cell.height

                        # LEFT SURFACE
                    
                        # transport term on diagonal
                        self.A[cn*ng + e, cn*ng + e] += cell.surfaces[0].DDif[e] * cell.length * cell.height
                        
                        # transport terms on off diagonals
                        if x != 0:
                            self.A[cn*ng + e, (cn - 1)*ng + e] -= cell.surfaces[0].DDif[e] * cell.length * cell.height

                        # FRONT SURFACE
                    
                        # transport term on diagonal
                        self.A[cn*ng + e, cn*ng + e] += cell.surfaces[3].DDif[e] * cell.width * cell.height
                        
                        # transport terms on off diagonals
                        if y != cy - 1:
                            self.A[cn*ng + e, (cn + cx)*ng + e] -= cell.surfaces[3].DDif[e] * cell.width * cell.height

                        # BACK SURFACE
                    
                        # transport term on diagonal
                        self.A[cn*ng + e, cn*ng + e] += cell.surfaces[2].DDif[e] * cell.width * cell.height
                        
                        # transport terms on off diagonals
                        if y != 0:
                            self.A[cn*ng + e, (cn - cx)*ng + e] -= cell.surfaces[2].DDif[e] * cell.width * cell.height

                        # BOTTOM SURFACE
                    
                        # transport term on diagonal
                        self.A[cn*ng + e, cn*ng + e] += cell.surfaces[4].DDif[e] * cell.width * cell.length
                        
                        # transport terms on off diagonals
                        if z != cz - 1:
                            self.A[cn*ng + e, (cn + cx*cy)*ng + e] -= cell.surfaces[4].DDif[e] * cell.width * cell.length

                        # TOP SURFACE
                    
                        # transport term on diagonal
                        self.A[cn*ng + e, cn*ng + e] += cell.surfaces[5].DDif[e] * cell.width * cell.length
                        
                        # transport terms on off diagonals
                        if z != 0:
                            self.A[cn*ng + e, (cn - cx*cy)*ng + e] -= cell.surfaces[5].DDif[e] * cell.width * cell.length

        self.A = self.A.tocsr()
        self.M = self.M.tocsr()
        
        
    def solve(self, tol, iterations):

        self.mesh.setCellBoundaries()
        self.computeDs()
        self.makeAM()
        self.computeFlux(tol)
        print 'DIFFUSION: --- k_eff = ' + str(self.keff)[0:10]


    def computeFlux(self, tol):
                 
        max_iter = 1000
        ng = self.mesh.cells[0].material.num_groups
        cx = self.mesh.cells_x
        cy = self.mesh.cells_y
        cz = self.mesh.cells_z
        print 'num cells x: ' + str(cx) + ', cy: ' + str(cy) + ', cz: ' + str(cz)
        sold = np.zeros(cx*cy*cz*ng)
        snew = np.zeros(cx*cy*cz*ng)
        res = np.zeros(cx*cy*cz*ng)
        self.phi = np.ones(cx*cy*cz*ng)
        self.keff_old = self.keff

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
        
        for i in range(max_iter):
            self.phi = spsolve(self.A, sold)
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
         
            # compute error
            l2_norm = l2_norm / (cx*cy*cz*ng)
            snew = snew * (cx*cy*cz*ng) / sumnew
             
            sold = snew
            
            print 'iteration: ' + str(i) + '  --- keff: ' + str(self.keff)[0:10]
             
            if l2_norm < tol:
                
                for i in range(cx*cy*cz):
                    for e in range(ng):
                        self.mesh.cells[i].flux[e] = self.phi[i*ng+e]
                break            
