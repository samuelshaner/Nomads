
import numpy as np
import matplotlib.pyplot as plt
from math import *
from Cell import *
from Material import *
from Mesh import *
import plotter as pttr
import os
import sys
import getopt
import time
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve


class Solver(object):

    def __init__(self, mesh, method):
        
        self.mesh = mesh
        self.method = method
        self.ng = mesh.cells[0].material.num_groups
        self.A = lil_matrix((mesh.num_cells*self.ng,mesh.num_cells*self.ng))
        self.M = lil_matrix((mesh.num_cells*self.ng,mesh.num_cells*self.ng))
        
        if method == 'NEM4':
            self.N = lil_matrix((8*mesh.num_cells*self.ng,8*mesh.num_cells*self.ng))
            self.coeffs = np.zeros(8*mesh.num_cells*self.ng)
            self.b = np.zeros(8*mesh.num_cells*self.ng)
        else:
            self.N = lil_matrix((4*mesh.num_cells*self.ng,4*mesh.num_cells*self.ng))
            self.coeffs = np.zeros(4*mesh.num_cells*self.ng)
            self.b = np.zeros(4*mesh.num_cells*self.ng)
            
        self.phi = np.ones(mesh.num_cells*self.ng)
        self.keff = 1.0
        self.keff_old = 0.0
        
    def computeDs(self):
        
        cw = self.mesh.cells_x
        ch = self.mesh.cells_y
        ng = self.ng
        flux_next = 0.0
        flux = 0.0
        
        relax_factor = 0.33
        
        for y in range(ch):
            for x in range(cw):
                
                cell = self.mesh.cells[y*cw+x]
                
                for side in range(4):
                    
                    surface = cell.surfaces[side]
                    cellNext = cell.neighborCells[side]
                    
                    for e in range(ng):
                        
                        d = cell.material.D[e]
                        flux = cell.flux[e]
                        
                        # set the sense of the surface
                        if side == 0 or side == 3:
                            sense = -1.0
                        else:
                            sense = 1.0
                            
                        # set the length of the surface parallel to and perpendicular from surface
                        if side == 0 or side == 2:
                            length_perpen = cell.width
                        elif side == 1 or side == 3:
                            length_perpen = cell.height
                            
                            
                        if cellNext is None:
                            if surface.boundary == 'reflective':
                                d_hat = 0.0
                                d_tilde = 0.0
                                current = 0.0
                                d_dif = 0.0
                                
                            else:
                                current = sense * surface.current[e]
                                d_dif = 2 * d / length_perpen / (1 + 4 * d / length_perpen)
                                d_hat = 2 * d / length_perpen / (1 + 4 * d / length_perpen)
                                if flux == 0.0:
                                    d_tilde = 0.0
                                else:
                                    d_tilde = (sense * d_hat * flux - current) / flux
                                
                        else:

                            if side == 0 or side == 2:
                                next_length_perpen = cellNext.width
                            elif side == 1 or side == 3:
                                next_length_perpen = cellNext.height

                            d_next = cellNext.material.D[e]
                            flux_next = cellNext.flux[e]
                        
                            d_dif = 2 * d * d_next / (length_perpen * d + next_length_perpen * d_next)
                            d_hat = 2 * d * d_next / (length_perpen * d + next_length_perpen * d_next)
                        
                            current = surface.current[e]
                            
                            if flux == 0.0 and flux_next == 0.0:
                                d_tilde = 0.0
                            else:
                                d_tilde = - (sense * d_hat * (flux_next - flux) + current) / (flux_next + flux)
                        
                        if abs(d_tilde) > abs(d_hat):
                            if sense == -1:
                                if (1 - abs(d_tilde)/d_tilde < 1.e-8):
                                    d_hat   = - current / (2*flux)
                                    d_tilde = - current / (2*flux)
                                else:
                                    d_hat   = current / (2*flux_next)
                                    d_tilde = - current / (2*flux_next)
                            else: 
                                if (1 - abs(d_tilde)/d_tilde < 1.e-8):
                                    d_hat   = - current / (2*flux_next)
                                    d_tilde = - current / (2*flux_next)
                                else:
                                    d_hat   = current / (2*flux)
                                    d_tilde = - current / (2*flux)
                                            
                                            
                        surface.DDif[e] = d_dif
                        surface.DHat[e] = d_hat
                        surface.DTilde[e] = surface.DTilde[e] * (1 - relax_factor) + relax_factor * d_tilde                
                         
                         
    def makeAM(self):
        
        cw = self.mesh.cells_x
        ch = self.mesh.cells_y
        ng = self.ng

        self.A = self.A.tolil()
        self.A = self.A * 0.0
        self.M = self.M.tolil()
        self.M = self.M * 0.0
        
        for y in range(ch): 
            for x in range(cw):
                
                for e in range(ng):
                    
                    cell = self.mesh.cells[y*cw+x]
                    
                    # absorption term on diagonal
                    self.A[(y*cw+x)*ng + e, (y*cw+x)*ng + e] += cell.material.sigma_a[e] * cell.volume
                    
                    # out scattering term on diagonal
                    for g in range(ng):
                        if e != g:
                            self.A[(y*cw+x)*ng + e, (y*cw+x)*ng + e] += cell.material.sigma_s[e,g] * cell.volume
                    
                    # fission terms on diagonal
                    for g in range(ng):
                        self.M[(y*cw+x)*ng + e, (y*cw+x)*ng + g] = cell.material.chi[e] * cell.material.nu_sigma_f[g] * cell.volume

                    # in scattering terms on off diagonals
                    for g in range(ng):
                        if e != g:
                            self.A[(y*cw+x)*ng + e, (y*cw+x)*ng + g] -= cell.material.sigma_s[g,e] * cell.volume
                            
                    # RIGHT SURFACE
                    
                    # transport term on diagonal
                    if self.method == 'NEM4' or self.method == 'NEM2':
                        self.A[(y*cw+x)*ng + e, (y*cw+x)*ng + e] += (cell.surfaces[2].DHat[e] - cell.surfaces[2].DTilde[e]) * cell.height
                    elif self.method == 'diffusion':
                        self.A[(y*cw+x)*ng + e, (y*cw+x)*ng + e] += cell.surfaces[2].DDif[e] * cell.height
                        
                    # transport terms on off diagonals
                    if x != cw - 1:
                        if self.method == 'NEM4' or self.method == 'NEM2':
                            self.A[(y*cw+x)*ng + e, (y*cw+x + 1)*ng + e] -= (cell.surfaces[2].DHat[e] + cell.surfaces[2].DTilde[e]) * cell.height
                        elif self.method == 'diffusion':
                            self.A[(y*cw+x)*ng + e, (y*cw+x + 1)*ng + e] -= cell.surfaces[2].DDif[e] * cell.height

                    # LEFT SURFACE
                    
                    # transport term on diagonal
                    if self.method == 'NEM4' or self.method == 'NEM2':
                        self.A[(y*cw+x)*ng + e, (y*cw+x)*ng + e] += (cell.surfaces[0].DHat[e] + cell.surfaces[0].DTilde[e]) * cell.height
                    elif self.method == 'diffusion':
                        self.A[(y*cw+x)*ng + e, (y*cw+x)*ng + e] += cell.surfaces[0].DDif[e] * cell.height
                        
                    # transport terms on off diagonals
                    if x != 0:
                        if self.method == 'NEM4' or self.method == 'NEM2':
                            self.A[(y*cw+x)*ng + e, (y*cw+x - 1)*ng + e] -= (cell.surfaces[0].DHat[e] - cell.surfaces[0].DTilde[e]) * cell.height
                        elif self.method == 'diffusion':
                            self.A[(y*cw+x)*ng + e, (y*cw+x - 1)*ng + e] -= cell.surfaces[0].DDif[e] * cell.height

                    # BOTTOM SURFACE
                    
                    # transport term on diagonal
                    if self.method == 'NEM4' or self.method == 'NEM2':
                        self.A[(y*cw+x)*ng + e, (y*cw+x)*ng + e] += (cell.surfaces[1].DHat[e] - cell.surfaces[1].DTilde[e]) * cell.width
                    elif self.method == 'diffusion':
                        self.A[(y*cw+x)*ng + e, (y*cw+x)*ng + e] += cell.surfaces[1].DDif[e] * cell.width
                        
                    # transport terms on off diagonals
                    if y != ch - 1:
                        if self.method == 'NEM4' or self.method == 'NEM2':
                            self.A[(y*cw+x)*ng + e, ((y+1)*cw+x)*ng + e] -= (cell.surfaces[1].DHat[e] + cell.surfaces[1].DTilde[e]) * cell.width
                        elif self.method == 'diffusion':
                            self.A[(y*cw+x)*ng + e, ((y+1)*cw+x)*ng + e] -= cell.surfaces[1].DDif[e] * cell.width

                    # TOP SURFACE
                    
                    # transport term on diagonal
                    if self.method == 'NEM4' or self.method == 'NEM2':
                        self.A[(y*cw+x)*ng + e, (y*cw+x)*ng + e] += (cell.surfaces[3].DHat[e] + cell.surfaces[3].DTilde[e]) * cell.width
                    elif self.method == 'diffusion':
                        self.A[(y*cw+x)*ng + e, (y*cw+x)*ng + e] += cell.surfaces[3].DDif[e] * cell.width
                        
                    # transport terms on off diagonals
                    if y != 0:
                        if self.method == 'NEM4' or self.method == 'NEM2':
                            self.A[(y*cw+x)*ng + e, ((y-1)*cw+x)*ng + e] -= (cell.surfaces[3].DHat[e] - cell.surfaces[3].DTilde[e]) * cell.width
                        elif self.method == 'diffusion':
                            self.A[(y*cw+x)*ng + e, ((y-1)*cw+x)*ng + e] -= cell.surfaces[3].DDif[e] * cell.width

        self.A = self.A.tocsr()
        self.M = self.M.tocsr()
                
           
    def makeN(self):
        
        cw = self.mesh.cells_x
        ch = self.mesh.cells_y
        ng = self.ng
        
        if self.method == 'NEM4':
            order = 4
        elif self.method == 'NEM2':
            order = 2
        
        self.N = self.N.tolil()
        self.N = self.N * 0.0
        
        for y in range(ch): 
            for x in range(cw):
                
                cell = self.mesh.cells[y*cw+x]
                
                for e in range(ng):
                    
                    cn = 2*order*e + 2*order*ng*(y*cw+x)
                    cr = 2*order*e + 2*order*ng*(y*cw+x+1)
                    cl = 2*order*e + 2*order*ng*(y*cw+x-1)
                    cd = 2*order*e + 2*order*ng*((y+1)*cw+x)
                    cu = 2*order*e + 2*order*ng*((y-1)*cw+x)                    
                    
                    if y == 0:
                        if x == 0:
                            # CURRENT BOUNDARY SURFACE 0
                            if order == 4:
                                self.N[cn,cn]   = self.dP1(0.0)
                                self.N[cn,cn+1] = self.dP2(0.0)
                                self.N[cn,cn+2] = self.dP3(0.0)
                                self.N[cn,cn+3] = self.dP4(0.0)
                            else:
                                self.N[cn,cn]   = self.dP1(0.0)
                                self.N[cn,cn+1] = self.dP2(0.0)
                            
                            if ch > 1:                             
                                cell_next = self.mesh.cells[(y+1)*cw+x]

                                # CURRENT INTERFACE SURFACE 1                        
                                if order == 4:
                                    self.N[cn+1,cn+4] = - self.dP1(1.0) * cell.material.D[e] / cell.height
                                    self.N[cn+1,cn+5] = - self.dP2(1.0) * cell.material.D[e] / cell.height
                                    self.N[cn+1,cn+6] = - self.dP3(1.0) * cell.material.D[e] / cell.height
                                    self.N[cn+1,cn+7] = - self.dP4(1.0) * cell.material.D[e] / cell.height
                                    self.N[cn+1,cd+4] = self.dP1(0.0) * cell_next.material.D[e] / cell_next.height
                                    self.N[cn+1,cd+5] = self.dP2(0.0) * cell_next.material.D[e] / cell_next.height
                                    self.N[cn+1,cd+6] = self.dP3(0.0) * cell_next.material.D[e] / cell_next.height
                                    self.N[cn+1,cd+7] = self.dP4(0.0) * cell_next.material.D[e] / cell_next.height                            
                                else:
                                    self.N[cn+1,cn+2] = - self.dP1(1.0) * cell.material.D[e] / cell.height
                                    self.N[cn+1,cn+3] = - self.dP2(1.0) * cell.material.D[e] / cell.height
                                    self.N[cn+1,cd+2] = self.dP1(0.0) * cell_next.material.D[e] / cell_next.height
                                    self.N[cn+1,cd+3] = self.dP2(0.0) * cell_next.material.D[e] / cell_next.height
                            else:
                                
                                # CURRENT BOUNDARY SURFACE 1
                                if order == 4:
                                    self.N[cn+1,cn+4] = self.dP1(1.0)
                                    self.N[cn+1,cn+5] = self.dP2(1.0)
                                    self.N[cn+1,cn+6] = self.dP3(1.0)
                                    self.N[cn+1,cn+7] = self.dP4(1.0)
                                else:
                                    self.N[cn+1,cn+2] = self.dP1(1.0)
                                    self.N[cn+1,cn+3] = self.dP2(1.0)
                       
                            # CURRENT INTERFACE SURFACE 2                        
                            if cw > 1:
                                cell_next = self.mesh.cells[y*cw+x+1]
                                
                                if order == 4:
                                    self.N[cn+2,cn]   = - self.dP1(1.0) * cell.material.D[e] / cell.width
                                    self.N[cn+2,cn+1] = - self.dP2(1.0) * cell.material.D[e] / cell.width
                                    self.N[cn+2,cn+2] = - self.dP3(1.0) * cell.material.D[e] / cell.width
                                    self.N[cn+2,cn+3] = - self.dP4(1.0) * cell.material.D[e] / cell.width
                                    self.N[cn+2,cr]   = self.dP1(0.0) * cell_next.material.D[e] / cell_next.width
                                    self.N[cn+2,cr+1] = self.dP2(0.0) * cell_next.material.D[e] / cell_next.width
                                    self.N[cn+2,cr+2] = self.dP3(0.0) * cell_next.material.D[e] / cell_next.width
                                    self.N[cn+2,cr+3] = self.dP4(0.0) * cell_next.material.D[e] / cell_next.width                            
                                else:
                                    self.N[cn+2,cn]   = - self.dP1(1.0) * cell.material.D[e] / cell.width
                                    self.N[cn+2,cn+1] = - self.dP2(1.0) * cell.material.D[e] / cell.width
                                    self.N[cn+2,cr]   = self.dP1(0.0) * cell_next.material.D[e] / cell_next.width
                                    self.N[cn+2,cr+1] = self.dP2(0.0) * cell_next.material.D[e] / cell_next.width
                            else:
                                # CURRENT BOUNDARY SURFACE 2
                                if order == 4:
                                    self.N[cn+2,cn]   = self.dP1(1.0)
                                    self.N[cn+2,cn+1] = self.dP2(1.0)
                                    self.N[cn+2,cn+2] = self.dP3(1.0)
                                    self.N[cn+2,cn+3] = self.dP4(1.0)
                                else:
                                    self.N[cn+2,cn]   = self.dP1(1.0)
                                    self.N[cn+2,cn+1] = self.dP2(1.0)
                        
                            # CURRENT BOUNDARY 3
                            if order == 4:
                                self.N[cn+3,cn+4] = self.dP1(0.0)
                                self.N[cn+3,cn+5] = self.dP2(0.0)
                                self.N[cn+3,cn+6] = self.dP3(0.0)
                                self.N[cn+3,cn+7] = self.dP4(0.0)
                            else:
                                self.N[cn+3,cn+2] = self.dP1(0.0)
                                self.N[cn+3,cn+3] = self.dP2(0.0)
                        
                        else:
                            
                            # FLUX SURFACE 0
                            self.N[cn,cn]   = - self.P1(0.0)
                            self.N[cn,cn+1] = - self.P2(0.0)
                            self.N[cn,cl]   = self.P1(1.0)
                            self.N[cn,cl+1] = self.P2(1.0)    
                            
                            self.b[cn] = - self.mesh.cells[y*cw+x-1].flux[e] + cell.flux[e]
                            
                            if ch > 1:                             
                                cell_next = self.mesh.cells[(y+1)*cw+x]

                                # CURRENT INTERFACE SURFACE 1                        
                                if order == 4:
                                    self.N[cn+1,cn+4] = - self.dP1(1.0) * cell.material.D[e] / cell.height
                                    self.N[cn+1,cn+5] = - self.dP2(1.0) * cell.material.D[e] / cell.height
                                    self.N[cn+1,cn+6] = - self.dP3(1.0) * cell.material.D[e] / cell.height
                                    self.N[cn+1,cn+7] = - self.dP4(1.0) * cell.material.D[e] / cell.height
                                    self.N[cn+1,cd+4] = self.dP1(0.0) * cell_next.material.D[e] / cell_next.height
                                    self.N[cn+1,cd+5] = self.dP2(0.0) * cell_next.material.D[e] / cell_next.height
                                    self.N[cn+1,cd+6] = self.dP3(0.0) * cell_next.material.D[e] / cell_next.height
                                    self.N[cn+1,cd+7] = self.dP4(0.0) * cell_next.material.D[e] / cell_next.height                            
                                else:
                                    self.N[cn+1,cn+2] = - self.dP1(1.0) * cell.material.D[e] / cell.height
                                    self.N[cn+1,cn+3] = - self.dP2(1.0) * cell.material.D[e] / cell.height
                                    self.N[cn+1,cd+2] = self.dP1(0.0) * cell_next.material.D[e] / cell_next.height
                                    self.N[cn+1,cd+3] = self.dP2(0.0) * cell_next.material.D[e] / cell_next.height
                            else:
                                
                                # CURRENT BOUNDARY SURFACE 1
                                if order == 4:
                                    self.N[cn+1,cn+4] = self.dP1(1.0)
                                    self.N[cn+1,cn+5] = self.dP2(1.0)
                                    self.N[cn+1,cn+6] = self.dP3(1.0)
                                    self.N[cn+1,cn+7] = self.dP4(1.0)
                                else:
                                    self.N[cn+1,cn+2] = self.dP1(1.0)
                                    self.N[cn+1,cn+3] = self.dP2(1.0)
                       
                            # CURRENT INTERFACE SURFACE 2                        
                            if x != cw - 1:
                                cell_next = self.mesh.cells[y*cw+x+1]
                                
                                if order == 4:
                                    self.N[cn+2,cn]   = - self.dP1(1.0) * cell.material.D[e] / cell.width
                                    self.N[cn+2,cn+1] = - self.dP2(1.0) * cell.material.D[e] / cell.width
                                    self.N[cn+2,cn+2] = - self.dP3(1.0) * cell.material.D[e] / cell.width
                                    self.N[cn+2,cn+3] = - self.dP4(1.0) * cell.material.D[e] / cell.width
                                    self.N[cn+2,cr]   = self.dP1(0.0) * cell_next.material.D[e] / cell_next.width
                                    self.N[cn+2,cr+1] = self.dP2(0.0) * cell_next.material.D[e] / cell_next.width
                                    self.N[cn+2,cr+2] = self.dP3(0.0) * cell_next.material.D[e] / cell_next.width
                                    self.N[cn+2,cr+3] = self.dP4(0.0) * cell_next.material.D[e] / cell_next.width                            
                                else:
                                    self.N[cn+2,cn]   = - self.dP1(1.0) * cell.material.D[e] / cell.width
                                    self.N[cn+2,cn+1] = - self.dP2(1.0) * cell.material.D[e] / cell.width
                                    self.N[cn+2,cr]   = self.dP1(0.0) * cell_next.material.D[e] / cell_next.width
                                    self.N[cn+2,cr+1] = self.dP2(0.0) * cell_next.material.D[e] / cell_next.width
                            else:
                                # CURRENT BOUNDARY SURFACE 2
                                if order == 4:
                                    self.N[cn+2,cn]   = self.dP1(1.0)
                                    self.N[cn+2,cn+1] = self.dP2(1.0)
                                    self.N[cn+2,cn+2] = self.dP3(1.0)
                                    self.N[cn+2,cn+3] = self.dP4(1.0)
                                else:
                                    self.N[cn+2,cn]   = self.dP1(1.0)
                                    self.N[cn+2,cn+1] = self.dP2(1.0)
                            
                            # CURRENT BOUNDARY 3
                            if order == 4:
                                self.N[cn+3,cn+4] = self.dP1(0.0)
                                self.N[cn+3,cn+5] = self.dP2(0.0)
                                self.N[cn+3,cn+6] = self.dP3(0.0)
                                self.N[cn+3,cn+7] = self.dP4(0.0)
                            else:
                                self.N[cn+3,cn+2] = self.dP1(0.0)
                                self.N[cn+3,cn+3] = self.dP2(0.0)
                                
                    elif y < cw - 1:
                        if x == 0:
                            # CURRENT BOUNDARY SURFACE 0
                            if order == 4:
                                self.N[cn,cn]   = self.dP1(0.0)
                                self.N[cn,cn+1] = self.dP2(0.0)
                                self.N[cn,cn+2] = self.dP3(0.0)
                                self.N[cn,cn+3] = self.dP4(0.0)
                            else:
                                self.N[cn,cn]   = self.dP1(0.0)
                                self.N[cn,cn+1] = self.dP2(0.0)
                        else: 
                            # FLUX SURFACE 0
                            self.N[cn,cn]   = - self.P1(0.0)
                            self.N[cn,cn+1] = - self.P2(0.0)
                            self.N[cn,cl]   = self.P1(1.0)
                            self.N[cn,cl+1] = self.P2(1.0)    
                            
                            self.b[cn] = - self.mesh.cells[y*cw+x-1].flux[e] + cell.flux[e]
                            
                        # CURRENT INTERFACE 1    
                        cell_next = self.mesh.cells[(y+1)*cw+x]

                        if order == 4:
                            self.N[cn+1,cn+4] = - self.dP1(1.0) * cell.material.D[e] / cell.height
                            self.N[cn+1,cn+5] = - self.dP2(1.0) * cell.material.D[e] / cell.height
                            self.N[cn+1,cn+6] = - self.dP3(1.0) * cell.material.D[e] / cell.height
                            self.N[cn+1,cn+7] = - self.dP4(1.0) * cell.material.D[e] / cell.height
                            self.N[cn+1,cd+4] = self.dP1(0.0) * cell_next.material.D[e] / cell_next.height
                            self.N[cn+1,cd+5] = self.dP2(0.0) * cell_next.material.D[e] / cell_next.height
                            self.N[cn+1,cd+6] = self.dP3(0.0) * cell_next.material.D[e] / cell_next.height
                            self.N[cn+1,cd+7] = self.dP4(0.0) * cell_next.material.D[e] / cell_next.height                            
                        else:
                            self.N[cn+1,cn+2] = - self.dP1(1.0) * cell.material.D[e] / cell.height
                            self.N[cn+1,cn+3] = - self.dP2(1.0) * cell.material.D[e] / cell.height
                            self.N[cn+1,cd+2] = self.dP1(0.0) * cell_next.material.D[e] / cell_next.height
                            self.N[cn+1,cd+3] = self.dP2(0.0) * cell_next.material.D[e] / cell_next.height
                        
                                                
                        if x != cw - 1:
                            cell_next = self.mesh.cells[y*cw+x+1]
                            
                            # CURRENT INTERFACE SURFACE 2
                            if order == 4:
                                self.N[cn+2,cn]   = - self.dP1(1.0) * cell.material.D[e] / cell.width
                                self.N[cn+2,cn+1] = - self.dP2(1.0) * cell.material.D[e] / cell.width
                                self.N[cn+2,cn+2] = - self.dP3(1.0) * cell.material.D[e] / cell.width
                                self.N[cn+2,cn+3] = - self.dP4(1.0) * cell.material.D[e] / cell.width
                                self.N[cn+2,cr]   = self.dP1(0.0) * cell_next.material.D[e] / cell_next.width
                                self.N[cn+2,cr+1] = self.dP2(0.0) * cell_next.material.D[e] / cell_next.width
                                self.N[cn+2,cr+2] = self.dP3(0.0) * cell_next.material.D[e] / cell_next.width
                                self.N[cn+2,cr+3] = self.dP4(0.0) * cell_next.material.D[e] / cell_next.width                            
                            else:
                                self.N[cn+2,cn]   = - self.dP1(1.0) * cell.material.D[e] / cell.width
                                self.N[cn+2,cn+1] = - self.dP2(1.0) * cell.material.D[e] / cell.width
                                self.N[cn+2,cr]   = self.dP1(0.0) * cell_next.material.D[e] / cell_next.width
                                self.N[cn+2,cr+1] = self.dP2(0.0) * cell_next.material.D[e] / cell_next.width
                        else:
                            # CURRENT BOUNDARY SURFACE 2
                            if order == 4:
                                self.N[cn+2,cn]   = self.dP1(1.0)
                                self.N[cn+2,cn+1] = self.dP2(1.0)
                                self.N[cn+2,cn+2] = self.dP3(1.0)
                                self.N[cn+2,cn+3] = self.dP4(1.0)
                            else:
                                self.N[cn+2,cn]   = self.dP1(1.0)
                                self.N[cn+2,cn+1] = self.dP2(1.0)                        
                    
                        # FLUX SURFACE 3
                        if order == 4:
                            self.N[cn+3,cn+4] = - self.P1(0.0)
                            self.N[cn+3,cn+5] = - self.P2(0.0)
                            self.N[cn+3,cu+4] = self.P1(1.0)
                            self.N[cn+3,cu+5] = self.P2(1.0)  
                        else:
                            self.N[cn+3,cn+2] = - self.P1(0.0)
                            self.N[cn+3,cn+3] = - self.P2(0.0)
                            self.N[cn+3,cu+2] = self.P1(1.0)
                            self.N[cn+3,cu+3] = self.P2(1.0)  
                        
                        self.b[cn+3] = - self.mesh.cells[(y-1)*cw+x].flux[e] + cell.flux[e]
                    
                    else:    
                        if x == 0:
                            # CURRENT BOUNDARY SURFACE 0
                            if order == 4:
                                self.N[cn,cn]   = self.dP1(0.0)
                                self.N[cn,cn+1] = self.dP2(0.0)
                                self.N[cn,cn+2] = self.dP3(0.0)
                                self.N[cn,cn+3] = self.dP4(0.0)
                            else:
                                self.N[cn,cn]   = self.dP1(0.0)
                                self.N[cn,cn+1] = self.dP2(0.0)
                        else: 
                            # FLUX SURFACE 0
                            self.N[cn,cn]   = - self.P1(0.0)
                            self.N[cn,cn+1] = - self.P2(0.0)
                            self.N[cn,cl]   = self.P1(1.0)
                            self.N[cn,cl+1] = self.P2(1.0)    
                            
                            self.b[cn] = - self.mesh.cells[y*cw+x-1].flux[e] + cell.flux[e]     
                                
                        # CURRENT BOUNDARY SURFACE 1
                        if order == 4:
                            self.N[cn+1,cn+4] = self.dP1(1.0)
                            self.N[cn+1,cn+5] = self.dP2(1.0)
                            self.N[cn+1,cn+6] = self.dP3(1.0)
                            self.N[cn+1,cn+7] = self.dP4(1.0)
                        else:
                            self.N[cn+1,cn+2] = self.dP1(1.0)
                            self.N[cn+1,cn+3] = self.dP2(1.0)   
                                
                        if x != cw - 1:
                            cell_next = self.mesh.cells[y*cw+x+1]
                            
                            # CURRENT INTERFACE SURFACE 2
                            if order == 4:
                                self.N[cn+2,cn]   = - self.dP1(1.0) * cell.material.D[e] / cell.width
                                self.N[cn+2,cn+1] = - self.dP2(1.0) * cell.material.D[e] / cell.width
                                self.N[cn+2,cn+2] = - self.dP3(1.0) * cell.material.D[e] / cell.width
                                self.N[cn+2,cn+3] = - self.dP4(1.0) * cell.material.D[e] / cell.width
                                self.N[cn+2,cr]   = self.dP1(0.0) * cell_next.material.D[e] / cell_next.width
                                self.N[cn+2,cr+1] = self.dP2(0.0) * cell_next.material.D[e] / cell_next.width
                                self.N[cn+2,cr+2] = self.dP3(0.0) * cell_next.material.D[e] / cell_next.width
                                self.N[cn+2,cr+3] = self.dP4(0.0) * cell_next.material.D[e] / cell_next.width                            
                            else:
                                self.N[cn+2,cn]   = - self.dP1(1.0) * cell.material.D[e] / cell.width
                                self.N[cn+2,cn+1] = - self.dP2(1.0) * cell.material.D[e] / cell.width
                                self.N[cn+2,cr]   = self.dP1(0.0) * cell_next.material.D[e] / cell_next.width
                                self.N[cn+2,cr+1] = self.dP2(0.0) * cell_next.material.D[e] / cell_next.width
                        else:
                            # CURRENT BOUNDARY SURFACE 2
                            if order == 4:
                                self.N[cn+2,cn]   = self.dP1(1.0)
                                self.N[cn+2,cn+1] = self.dP2(1.0)
                                self.N[cn+2,cn+2] = self.dP3(1.0)
                                self.N[cn+2,cn+3] = self.dP4(1.0)
                            else:
                                self.N[cn+2,cn]   = self.dP1(1.0)
                                self.N[cn+2,cn+1] = self.dP2(1.0)                                 
                                
                        # FLUX SURFACE 3
                        if order == 4:
                            self.N[cn+3,cn+4] = - self.P1(0.0)
                            self.N[cn+3,cn+5] = - self.P2(0.0)
                            self.N[cn+3,cu+4] = self.P1(1.0)
                            self.N[cn+3,cu+5] = self.P2(1.0)  
                        else:
                            self.N[cn+3,cn+2] = - self.P1(0.0)
                            self.N[cn+3,cn+3] = - self.P2(0.0)
                            self.N[cn+3,cu+2] = self.P1(1.0)
                            self.N[cn+3,cu+3] = self.P2(1.0)  
                        
                        self.b[cn+3] = - self.mesh.cells[(y-1)*cw+x].flux[e] + cell.flux[e]                                                            

                    if order == 4:
                        
                        # RESIDUAL X1
                        self.N[cn+4,cn]   = cell.material.sigma_r[e] / 3.0
                        self.N[cn+4,cn+1] = 0.0
                        self.N[cn+4,cn+2] = cell.material.sigma_r[e] / 5.0 + 12 * cell.material.D[e] / cell.width**2
                        self.N[cn+4,cn+3] = 0.0
                        
                        cnp = 2*order*ng*(y*cw+x)
                            
                        # in scattering and fission
                        for g in range(ng):
                            # fission
                            self.N[cn+4,cnp+2*order*g]   -= cell.material.chi[e] * cell.material.nu_sigma_f[g] / self.keff / 3.0
                            self.N[cn+4,cnp+2*order*g+2] -= cell.material.chi[e] * cell.material.nu_sigma_f[g] / self.keff / 5.0
                                
                            # in scattering
                            if g != e:
                                self.N[cn+4,cnp+2*order*g]   -= cell.material.sigma_s[g,e] / 3.0
                                self.N[cn+4,cnp+2*order*g+2] -= cell.material.sigma_s[g,e] / 5.0

                            
                        # RESIDUAL X2
                        self.N[cn+5,cn]   = 0.0
                        self.N[cn+5,cn+1] = cell.material.sigma_r[e] / 5.0
                        self.N[cn+5,cn+2] = 0.0
                        self.N[cn+5,cn+3] = - 3.0 * cell.material.sigma_r[e] / 35.0 - 12.0 * cell.material.D[e] / cell.width**2
                            
                        # in scattering and fission
                        for g in range(ng):
                            # fission
                            self.N[cn+5,cnp+2*order*g+1] -= cell.material.chi[e] * cell.material.nu_sigma_f[g] / self.keff / 5.0
                            self.N[cn+5,cnp+2*order*g+3] -= cell.material.chi[e] * cell.material.nu_sigma_f[g] / self.keff * (-3.0) / 35.0
                                
                            # in scattering
                            if g != e:
                                self.N[cn+5,cnp+2*order*g+1] -= cell.material.sigma_s[g,e] / 5.0
                                self.N[cn+5,cnp+2*order*g+3] -= cell.material.sigma_s[g,e] * (-3.0) / 35.0

                        # RESIDUAL Y1
                        self.N[cn+6,cn+4] = cell.material.sigma_r[e] / 3.0
                        self.N[cn+6,cn+5] = 0.0
                        self.N[cn+6,cn+6] = cell.material.sigma_r[e] / 5.0 + 12 * cell.material.D[e] / cell.height**2
                        self.N[cn+6,cn+7] = 0.0

                        # in scattering and fission
                        for g in range(ng):
                            # fission
                            self.N[cn+6,cnp+2*order*g+4] -= cell.material.chi[e] * cell.material.nu_sigma_f[g] / self.keff / 3.0
                            self.N[cn+6,cnp+2*order*g+6] -= cell.material.chi[e] * cell.material.nu_sigma_f[g] / self.keff / 5.0
                                
                            # in scattering
                            if g != e:
                                self.N[cn+6,cnp+2*order*g+4] -= cell.material.sigma_s[g,e] / 3.0
                                self.N[cn+6,cnp+2*order*g+6] -= cell.material.sigma_s[g,e] / 5.0
                          
                        # RESIDUAL Y2
                        self.N[cn+7,cn+4] = 0.0
                        self.N[cn+7,cn+5] = cell.material.sigma_r[e] / 5.0
                        self.N[cn+7,cn+6] = 0.0
                        self.N[cn+7,cn+7] = - 3.0 * cell.material.sigma_r[e] / 35.0 - 12.0 * cell.material.D[e] / cell.height**2
                            
                        # in scattering and fission
                        for g in range(ng):
                            # fission
                            self.N[cn+7,cnp+2*order*g+5] -= cell.material.chi[e] * cell.material.nu_sigma_f[g] / self.keff / 5.0
                            self.N[cn+7,cnp+2*order*g+7] -= cell.material.chi[e] * cell.material.nu_sigma_f[g] / self.keff * (-3.0) / 35.0
                                
                            # in scattering
                            if g != e:
                                self.N[cn+7,cnp+2*order*g+5] -= cell.material.sigma_s[g,e] / 5.0
                                self.N[cn+7,cnp+2*order*g+7] -= cell.material.sigma_s[g,e] * (-3.0) / 35.0
                        
                        
        self.N = self.N.tocsr()

        
    def computeCoeffs(self):
        
        self.coeffs = spsolve(self.N,self.b)
        self.N = self.N.tolil()
        
    # compute the currents        
    def computeCurrents(self):
        
        cw = self.mesh.cells_x
        ch = self.mesh.cells_y
        ng = self.ng
        
        if self.method == 'NEM4':
            order = 4
        elif self.method == 'NEM2':
            order = 2
        
        for y in range(ch): 
            for x in range(cw):
                
                cell = self.mesh.cells[y*cw+x]
                
                for e in range(ng):
                    
                    for side in range(4):
                        surface = cell.surfaces[side]
                        
                        if cell.neighborCells[side] is not None:
                            cn = 2*order*ng*(y*cw+x) + 2*order*e
                            
                            if order == 4:
                                if side == 0:
                                    surface.current[e] = - cell.material.D[e] / cell.width  * (self.coeffs[cn]   * self.dP1(0.0) + self.coeffs[cn+1] * self.dP2(0.0) + self.coeffs[cn+2] * self.dP3(0.0) + self.coeffs[cn+3] * self.dP4(0.0))
                                elif side == 1:
                                    surface.current[e] = - cell.material.D[e] / cell.height * (self.coeffs[cn+4] * self.dP1(1.0) + self.coeffs[cn+5] * self.dP2(1.0) + self.coeffs[cn+6] * self.dP3(1.0) + self.coeffs[cn+7] * self.dP4(1.0))
                                elif side == 2:
                                    surface.current[e] = - cell.material.D[e] / cell.width  * (self.coeffs[cn]   * self.dP1(1.0) + self.coeffs[cn+1] * self.dP2(1.0) + self.coeffs[cn+2] * self.dP3(1.0) + self.coeffs[cn+3] * self.dP4(1.0))
                                else:
                                    surface.current[e] = - cell.material.D[e] / cell.height * (self.coeffs[cn+4] * self.dP1(0.0) + self.coeffs[cn+5] * self.dP2(0.0) + self.coeffs[cn+6] * self.dP3(0.0) + self.coeffs[cn+7] * self.dP4(0.0))
                                                                        
                            elif order == 2:
                                if side == 0:
                                    surface.current[e] = - cell.material.D[e] / cell.width  * (self.coeffs[cn]   * self.dP1(0.0) + self.coeffs[cn+1] * self.dP2(0.0))
                                elif side == 1:
                                    surface.current[e] = - cell.material.D[e] / cell.height * (self.coeffs[cn+2] * self.dP1(1.0) + self.coeffs[cn+3] * self.dP2(1.0))
                                elif side == 2:
                                    surface.current[e] = - cell.material.D[e] / cell.width * (self.coeffs[cn]    * self.dP1(1.0) + self.coeffs[cn+1] * self.dP2(1.0))
                                else:
                                    surface.current[e] = - cell.material.D[e] / cell.height * (self.coeffs[cn+2] * self.dP1(0.0) + self.coeffs[cn+3] * self.dP2(0.0))
                                    
    def P1(self, val):
        return (2 * val - 1)
    
    def dP1(self, val):
        return 2
    
    def P2(self, val):
        return (6.0*val*(1-val)-1)
    
    def dP2(self, val):
        return (6 - 12.0*val)
    
    def P3(self, val):
        return (6.0*val*(1-val)*(2*val-1))

    def dP3(self, val):
        return (-6.0*(1-6*val+6*val**2))
    
    def P4(self, val):
        return (6.0*val*(1-val)*(5*val**2-5*val+1))
    
    def dP4(self, val):
        return (6.0-72*val+180*val**2-120*val**3)     
        

    def solve(self, tol):
                 
        max_iter = 1000
        ng = self.mesh.cells[0].material.num_groups
        cw = self.mesh.cells_x
        ch = self.mesh.cells_y
        sold = np.zeros(cw*ch*ng)
        snew = np.zeros(cw*ch*ng)
        res = np.zeros(cw*ch*ng)
        self.phi = np.ones(cw*ch*ng)
        self.keff_old = self.keff

        # get initial source and find initial keff
        snew = self.M * self.phi
        sumnew = np.sum(snew)
        sold = self.A * self.phi
        sumold = np.sum(sold)
        keff = sumnew / sumold
        
        # recompute and initialize the initial source
        sold = self.M * self.phi
        sumold = sum(sold)
        sold = sold * (cw*ch*ng / sumold)
        sumold = cw*ch*ng
        
        for i in range(max_iter):
            self.phi = spsolve(self.A, sold)
            snew = self.M * self.phi
            sumnew = np.sum(snew)
            
            # compute and set keff
            keff = sumnew / sumold
            sold = sold * keff
             
            # compute L2 norm of source error
            snew += 1.e-10
            res = np.divide(sold, snew)
            res = res - 1.0
            l2_norm = np.linalg.norm(res)
         
            # compute error
            l2_norm = l2_norm / (cw*ch*ng)
            snew = snew * (cw*ch*ng) / sumnew
             
            sold = snew
            
            print 'iteration: ' + str(i) + '  --- keff: ' + str(keff)[0:10]
             
            if l2_norm < tol:
                self.keff = keff
                
                for i in range(cw*ch):
                    for e in range(ng):
                        self.mesh.cells[i].flux[e] = self.phi[i*ng+e]
                break            


def main():
    
    start = time.time()
    
    # parse command line options
    try:
        opts, args = getopt.getopt(sys.argv[1:], "t:i:f:m:w:h:c:", ["tolerance","iterations","flux_plot", "method", "mesh_width", "mesh_height", "cell_size"])
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(2)

    # set default values
    tol = 1.e-8
    plot_flux = False
    
    mesh_width   = 2
    mesh_height  = 1
    cell_size    = 10.0
    fuel_width   = 10.0
    solve_method = 'NEM4'
    iterations = 50

    for o, a in opts:
        if o in ("-t", "--tolerance"):
            tol = float(a)
        elif o in ("-f", "--flux_plot"):
            plot_flux = True
        elif o in ("-m", "--method"):
            solve_method = a
        elif o in ("-w", "--mesh_width"):
            mesh_width = int(a)
        elif o in ("-h", "--mesh_height"):
            mesh_height = int(a)
        elif o in ("-c", "--cell_size"):
            cell_size = int(a)
        elif o in ("-i", "--iterations"):
            iterations = int(a)

        else:
            assert False, "unhandled option"

    # create mesh
    mesh = Mesh([cell_size,cell_size], [cell_size])
      
    # create fuel
    fuel = Material(2, 'fuel')
    fuel.setSigmaA([0.005, 0.10])
    fuel.setD([1.5, 0.40])
    fuel.setNuSigmaF([0.005, 0.15])
    fuel.setChi([1.0, 0.0])
    fuel.setSigmaS(np.array([[0.0, 0.02],[0.0, 0.0]]))
    
    # create fuel
    moderator = Material(2, 'moderator')
    moderator.setSigmaA([0.0, 0.01])
    moderator.setD([1.5, 0.20])
    moderator.setSigmaS(np.array([[0.0, 0.025],[0.0, 0.0]]))
    
    if solve_method == 'NEM4':
        order = 4
    else:
        order = 2
        
    mesh.cells[0].setMaterial(fuel, order)
    mesh.cells[1].setMaterial(moderator, order)
    mesh.makeSurfaces()
    
    # plot the mesh
    pttr.plotMesh(mesh)
    
    # create solver
    solver = Solver(mesh, solve_method)   

    # solve the matrix problem to get flux profile and keff
    
    if solve_method == 'NEM4' or solve_method == 'NEM2':
        for iteration in range(iterations):
        
            print 'CMFD outer iteration ' + str(iteration)
        
            solver.computeDs()
            solver.makeAM()
            solver.solve(tol)
            solver.makeN()
            solver.computeCoeffs()
            solver.computeCurrents()
        
            if abs(solver.keff_old - solver.keff) < 1.e-8:
                print solve_method + ': Converged in ' + str(iteration+1) + ' iterations --- k_eff = ' + str(solver.keff)[0:10]
                break
                                
    elif solve_method == 'diffusion':
        solver.computeDs()
        solver.makeAM()
        solver.solve(tol)
        print 'DIFFUSION: --- k_eff = ' + str(solver.keff)[0:10]


    
    if solve_method == 'NEM4' or solve_method == 'NEM2':
        pttr.plotFlux(solver)
        pttr.plotCellFlux(solver)
#         pttr.plotCurrent(solver)

    stop = time.time()
    
    print 'Ran time ' + str(stop-start)[0:5] + ' seconds'

    print '----------------------------------------------------------------------'

if __name__ == '__main__':

    main()



