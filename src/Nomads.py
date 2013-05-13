
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


class Solver(object):

    def __init__(self, mesh, method):
        
        self.mesh = mesh
        self.method = method
        self.ng = mesh.cells[0].material.num_groups
        self.A = np.zeros([mesh.num_cells*self.ng,mesh.num_cells*self.ng])
        self.M = np.zeros([mesh.num_cells*self.ng,mesh.num_cells*self.ng])
        self.phi = np.ones(mesh.num_cells*self.ng)
        
    def computeDs(self):
        
        cw = self.mesh.cells_x
        ch = self.mesh.cells_y
        ng = self.ng
        flux_next = 0.0
        flux = 0.0
        
        for y in range(ch):
            for x in range(cw):
                
                cell = self.mesh.cells[y*cw+x]
                
                for side in range(4):
                    
                    surface = cell.surfaces[side]
                    cellNext = cell.neighborCells[side]
                    
                    for e in range(ng):
                        
                        d = cell.material.D[e]
                        flux = cell.old_flux[e]
                        
                        # set the sense of the surface
                        if side == 0 or side == 3:
                            sense = -1.0
                        else:
                            sense = 1.0
                            
                        # set the length of the surface parallel to and perpendicular from surface
                        if side == 0 or side == 2:
                            length = cell.height
                            length_perpen = cell.width
                        elif side == 1 or side == 3:
                            length = cell.width
                            length_perpen = cell.height
                            
                        if cellNext == '':
                            if surface.boundary == 'reflective':
                                d_hat = 0.0
                                d_tilde = 0.0
                                current = 0.0
                                
                                surface.DDif[e] = 0.0
                                
                            else:
                                current = sense * surface.current[e]
                                surface.DDif[e] = 2 * d / length_perpen / (1 + 4 * d / length_perpen)
                                d_hat = 2 * d / length_perpen / (1 + 4 * d / length_perpen)
                                d_tilde = (sense * d_hat * flux - current / length) / flux
                                
                        else:
                            if side == 0:
                                next_length_perpen = cellNext.width
                                next_side = 2
                            elif side == 1:
                                next_length_perpen = cellNext.height
                                next_side = 3
                            elif side == 2:
                                next_length_perpen = cellNext.width
                                next_side = 0
                            elif side == 3:
                                next_length_perpen = cellNext.height
                                next_side = 1

                            d_next = cellNext.material.D[e]
                            flux_next = cellNext.old_flux[e]
                        
                            surface.DDif[e] = 2 * d * d_next / (length_perpen * d + next_length_perpen * d_next)
                        
                            d_hat = 2 * d * d_next / (length_perpen * d + next_length_perpen * d_next)
                        
                            current = 0.0
                        
                            current += sense * surface.current[e]
                        
                            current -= sense * cellNext.surfaces[next_side].current[e]
                        
                            d_tilde = - (sense * d_hat * (flux_next - flux) + current / length) / (flux_next + flux)
                        
                        if abs(d_tilde) > abs(d_hat):
                        
                            if sense == -1:
                                if (1 - abs(d_tilde)/d_tilde < 1.e-8):
                                    d_hat   = - current / (2*flux*length)
                                    d_tilde = - current / (2*flux*length)
                                else:
                                    d_hat   = current / (2*flux_next*length)
                                    d_tilde = - current / (2*flux_next*length)
                            else: 
                                if (1 - abs(d_tilde)/d_tilde < 1.e-8):
                                    d_hat   = - current / (2*flux_next*length)
                                    d_tilde = - current / (2*flux_next*length)
                                else:
                                    d_hat   = current / (2*flux*length)
                                    d_tilde = - current / (2*flux*length)
                                            
                                            
                        surface.DHat[e] = d_hat
                        surface.DTilde[e] = d_tilde
                
                
                
                
                
                
                         
    def makeAM(self):
        
        cw = self.mesh.cells_x
        ch = self.mesh.cells_y
        ng = self.ng
        
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
                    if self.method == 'cmfd':
                        self.A[(y*cw+x)*ng + e, (y*cw+x)*ng + e] += (cell.surfaces[2].DHat[e] - cell.surfaces[2].DTilde[e]) * cell.height
                    elif self.method == 'diffusion':
                        self.A[(y*cw+x)*ng + e, (y*cw+x)*ng + e] += cell.surfaces[2].DDif[e] * cell.height
                        
                    # transport terms on off diagonals
                    if x != cw - 1:
                        if self.method == 'cmfd':
                            self.A[(y*cw+x)*ng + e, (y*cw+x + 1)*ng + e] -= (cell.surfaces[2].DHat[e] + cell.surfaces[2].DTilde[e]) * cell.height
                        elif self.method == 'diffusion':
                            self.A[(y*cw+x)*ng + e, (y*cw+x + 1)*ng + e] -= cell.surfaces[2].DDif[e] * cell.height

                    # LEFT SURFACE
                    
                    # transport term on diagonal
                    if self.method == 'cmfd':
                        self.A[(y*cw+x)*ng + e, (y*cw+x)*ng + e] += (cell.surfaces[0].DHat[e] + cell.surfaces[0].DTilde[e]) * cell.height
                    elif self.method == 'diffusion':
                        self.A[(y*cw+x)*ng + e, (y*cw+x)*ng + e] += cell.surfaces[0].DDif[e] * cell.height
                        
                    # transport terms on off diagonals
                    if x != 0:
                        if self.method == 'cmfd':
                            self.A[(y*cw+x)*ng + e, (y*cw+x - 1)*ng + e] -= (cell.surfaces[0].DHat[e] - cell.surfaces[0].DTilde[e]) * cell.height
                        elif self.method == 'diffusion':
                            self.A[(y*cw+x)*ng + e, (y*cw+x - 1)*ng + e] -= cell.surfaces[0].DDif[e] * cell.height

                    # BOTTOM SURFACE
                    
                    # transport term on diagonal
                    if self.method == 'cmfd':
                        self.A[(y*cw+x)*ng + e, (y*cw+x)*ng + e] += (cell.surfaces[1].DHat[e] - cell.surfaces[1].DTilde[e]) * cell.width
                    elif self.method == 'diffusion':
                        self.A[(y*cw+x)*ng + e, (y*cw+x)*ng + e] += cell.surfaces[1].DDif[e] * cell.width
                        
                    # transport terms on off diagonals
                    if y != ch - 1:
                        if self.method == 'cmfd':
                            self.A[(y*cw+x)*ng + e, ((y+1)*cw+x)*ng + e] -= (cell.surfaces[1].DHat[e] + cell.surfaces[1].DTilde[e]) * cell.width
                        elif self.method == 'diffusion':
                            self.A[(y*cw+x)*ng + e, ((y+1)*cw+x)*ng + e] -= cell.surfaces[1].DDif[e] * cell.width

                    # TOP SURFACE
                    
                    # transport term on diagonal
                    if self.method == 'cmfd':
                        self.A[(y*cw+x)*ng + e, (y*cw+x)*ng + e] += (cell.surfaces[3].DHat[e] + cell.surfaces[3].DTilde[e]) * cell.width
                    elif self.method == 'diffusion':
                        self.A[(y*cw+x)*ng + e, (y*cw+x)*ng + e] += cell.surfaces[3].DDif[e] * cell.width
                        
                    # transport terms on off diagonals
                    if y != 0:
                        if self.method == 'cmfd':
                            self.A[(y*cw+x)*ng + e, ((y-1)*cw+x)*ng + e] -= (cell.surfaces[3].DHat[e] - cell.surfaces[3].DTilde[e]) * cell.width
                        elif self.method == 'diffusion':
                            self.A[(y*cw+x)*ng + e, ((y-1)*cw+x)*ng + e] -= cell.surfaces[3].DDif[e] * cell.width


        
       
#     # compute the NEM currents
#     def computeCurrents(self):
#          
        
        
    def solve(self, tol):
         
        max_iter = 10
        ng = self.mesh.cells[0].material.num_groups
        cw = self.mesh.cells_x
        ch = self.mesh.cells_y
        sold = np.zeros(cw*ch*ng)
        snew = np.zeros(cw*ch*ng)
        res = np.zeros(cw*ch*ng)
        phi_new = np.ones(cw*ch*ng)
        criteria = 1.e-5

        # get initial source and find initial keff
        snew = np.dot(self.M,phi_new)
        sumnew = np.sum(snew)
        sold = np.dot(self.A,phi_new)
        sumold = np.sum(sold)
        keff = sumnew / sumold
        
        # recompute and initialize the initial source
        sold = np.dot(self.M,phi_new)
        sumold = np.sum(sold)
        sold = sold * (cw*ch*ng / sumold)
        sumold = cw*ch*ng
        
        for i in range(max_iter):
            phi_new = np.linalg.solve(self.A, sold)
            snew = np.dot(self.M,phi_new)
            sumnew = np.sum(snew)
            
            # compute and set keff
            keff = sumnew / sumold
            sold = sold * keff
             
            # compute L2 norm of source error
            snew += 1.e-10
            res = np.divide(sold, snew)
            res = res - 1.0
            l2_norm = 0.0
            
            for j in res:
                l2_norm += j**2
                 
            l2_norm = sqrt(l2_norm) 
         
            # compute error
            l2_norm = l2_norm / (cw*ch*ng)
            snew = snew * (cw*ch*ng) / sumnew
             
            sold = snew
            
            print 'iteration: ' + str(i) + '  --- keff: ' + str(keff)
             
             
            if i > 5 and l2_norm < criteria:
                break            
            



def main():
    
    # parse command line options
    try:
        opts, args = getopt.getopt(sys.argv[1:], "t:f:m:w:h:c:", ["tolerance","flux_plot", "method", "mesh_width", "mesh_height", "cell_size"])
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(2)

    # set default values
    tol = .0001
    plot_flux = False
    
    mesh_width   = 2
    mesh_height  = 1
    cell_size    = 10.0
    fuel_width   = 10.0
    solve_method = 'diffusion'

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

        else:
            assert False, "unhandled option"

    # create mesh
    mesh = Mesh([cell_size,cell_size], [1.0])
       
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

    # assign fuel and moderator materials to mesh
    mesh.cells[0].setMaterial(fuel)
    mesh.cells[1].setMaterial(moderator)
    
    # plot the mesh
    pttr.plotMesh(mesh)
    
    # create solver
    solver = Solver(mesh, solve_method)   
    
    start = time.time()
    solver.computeDs()
    solver.makeAM()

    # solve the matrix problem to get flux profile and keff
    solver.solve(tol)

    stop = time.time()

    print 'Ran solver in ' + str(stop-start) + ' seconds'

#     if plot_flux:
#         pttr.plotFlux(mesh)

    print '----------------------------------------------------------------------'

if __name__ == '__main__':

    main()



