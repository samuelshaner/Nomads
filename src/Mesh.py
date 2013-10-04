
from math import *
import numpy as np
from Cell import *
from Surface import *
import sys


class Mesh(object):

    def __init__(self, lengths_x, lengths_y, lengths_z=[]):
        
        # set the properties of the mesh cells
        self.cells_x   = len(lengths_x)
        self.cells_y   = len(lengths_y)
        self.lengths_x = lengths_x
        self.lengths_y = lengths_y
        self.length_x  = sum(lengths_x)
        self.length_y  = sum(lengths_y)

        # create mesh boundaries
        self.boundaries = np.empty(6, dtype=object)
        for i in range(6):
            self.boundaries[i] = Surface(0)

        # check if user input info for 3rd dimension
        if lengths_z:
            self.dimensions = 3
            self.lengths_z  = lengths_z
            self.length_z   = sum(lengths_z)
            self.cells_z    = len(lengths_z)
        else:
            self.dimensions = 2
            self.lengths_z  = [1.0]
            self.length_z   = 1.0
            self.cells_z    = 1

        # create array of Cell object pointers and set number of cells
        self.num_cells = self.cells_x * self.cells_y * self.cells_z
        self.cells = np.empty(self.num_cells, dtype=object)

        # create Cell objects
        for z in range(self.cells_z):
            for y in range(self.cells_y):
                for x in range(self.cells_x):
                    self.cells[z*(self.cells_x*self.cells_y) + y*self.cells_x + x] = Cell(self.lengths_x[x], self.lengths_y[y], self.lengths_z[z])
                                
        # set neighbor cells
        for z in range(self.cells_z):
            for y in range(self.cells_y):
                for x in range(self.cells_x):

                    cell_num = z*(self.cells_x*self.cells_y) + y*self.cells_x + x
                    cell = self.cells[cell_num]

                    if x != 0:
                        cell.neighborCells[0] = self.cells[cell_num - 1]
         
                    if x != self.cells_x - 1:
                        cell.neighborCells[1] = self.cells[cell_num + 1]

                    if y != 0:
                        cell.neighborCells[2] = self.cells[cell_num - self.cells_x]
        
                    if y != self.cells_y - 1:
                        cell.neighborCells[3] = self.cells[cell_num + self.cells_x]

                    if z != 0:
                        cell.neighborCells[4] = self.cells[cell_num - self.cells_x*self.cells_y]
        
                    if z != self.cells_z - 1:
                        cell.neighborCells[5] = self.cells[cell_num + self.cells_x*self.cells_y]
                    

    # set the mesh boundary conditions (VACUUM, ZERO_FLUX, or REFLECTIVE)
    def setBoundaries(self, left, right, back, front, bottom='REFLECTIVE', top='REFLECTIVE'):

        self.boundaries[0].boundary = left
        self.boundaries[1].boundary = right
        self.boundaries[2].boundary = back
        self.boundaries[3].boundary = front
        self.boundaries[4].boundary = bottom
        self.boundaries[5].boundary = top

    # make mesh surfaces for each cell
    def makeSurfaces(self):

        for z in range(self.cells_z):
            for y in range(self.cells_y):
                for x in range(self.cells_x):
         
                    cell_num = z*self.cells_x*self.cells_y + y*self.cells_x + x
                    cell = self.cells[cell_num]

                    if x == 0:
                        cell.surfaces[0] = Surface(cell.material.num_groups)
                    else:
                        cell.surfaces[0] = self.cells[cell_num - 1].surfaces[1]

                    if y == 0:
                        cell.surfaces[2] = Surface(cell.material.num_groups)
                    else:
                        cell.surfaces[2] = self.cells[cell_num - self.cells_x].surfaces[3]

                    if z == 0:
                        cell.surfaces[4] = Surface(cell.material.num_groups)
                    else:
                        cell.surfaces[4] = self.cells[cell_num - self.cells_x*self.cells_y].surfaces[5]

                    cell.surfaces[1] = Surface(cell.material.num_groups)
                
                    cell.surfaces[3] = Surface(cell.material.num_groups)

                    cell.surfaces[5] = Surface(cell.material.num_groups)                
                    
    
    # break the geometry up into very small cells and solve diffusion problem
    def refineMesh(self, cell_size, dimensions=['x','y']):
        
        # check if mesh is being refined in the x direction
        if 'x' in dimensions:

            # break up the old lengths_x into chunks
            lengths_x = []
            for i in self.lengths_x:
                new_length_x = int(floor(i/cell_size))
                new_lengths_x = []

                # set new cell lengths
                for j in range(new_length_x):
                    new_lengths_x.append(cell_size)
                
                # if new cells don't entirely fill old cell,
                # add a small cell to fill the gap
                if abs(sum(new_lengths_x) - i) >= 1.e-3:
                    new_lengths_x.append(i - sum(new_lengths_x))
                
                for k in new_lengths_x:
                    lengths_x.append(k) 
        else:
            lengths_x = self.lengths_x

        if 'y' in dimensions:
 
           # break up the old lengths_y into chunks
            lengths_y = []
            for i in self.lengths_y:
                new_length_y = int(floor(i/cell_size))
                new_lengths_y = []
            
                for j in range(new_length_y):
                    new_lengths_y.append(cell_size)
                
                if abs(sum(new_lengths_y) - i) >= 1.e-3:
                    new_lengths_y.append(i - sum(new_lengths_y))
                
                for k in new_lengths_y:
                    lengths_y.append(k)

        else:
            lengths_y = self.lengths_y

        if 'z' in dimensions:

            # break up the old lengths_z into chunks
            lengths_z = []
            for i in self.lengths_z:
                new_length_z = int(floor(i/cell_size))
                new_lengths_z = []
            
                for j in range(new_length_z):
                    new_lengths_z.append(cell_size)
                
                if abs(sum(new_lengths_z) - i) >= 1.e-3:
                    new_lengths_z.append(i - sum(new_lengths_z))
                
                for k in new_lengths_z:
                    lengths_z.append(k)
        else:
            lengths_z = self.lengths_z

        if self.dimensions == 2:
            mesh = Mesh(lengths_x, lengths_y)
        else:
            mesh = Mesh(lengths_x, lengths_y, lengths_z)

        # set the cell materials
        for zn in range(mesh.cells_z):
            for yn in range(mesh.cells_y):
                for xn in range(mesh.cells_x):
            
                    # get center coordinate of new cell
                    xval = (sum(mesh.lengths_x[0:xn]) + sum(mesh.lengths_x[0:xn+1]))/2.0
                    yval = (sum(mesh.lengths_y[0:yn]) + sum(mesh.lengths_y[0:yn+1]))/2.0
                    zval = (sum(mesh.lengths_z[0:zn]) + sum(mesh.lengths_z[0:zn+1]))/2.0
           
                    # find which old cell the new cell is in
                    flag = 0     
                    for x in range(len(self.lengths_x)):
                        for y in range(len(self.lengths_y)):
                            for z in range(len(self.lengths_z)):
                                
                                cell_num = z*self.cells_x*self.cells_y + y*self.cells_x + x
                                cell = self.cells[cell_num]
                                
                                # if new cell center is bounded by old cell, set the new cell's material
                                if (xval >= sum(self.lengths_x[0:x]) and xval <= sum(self.lengths_x[0:x+1]) and 
                                    yval >= sum(self.lengths_y[0:y]) and yval <= sum(self.lengths_y[0:y+1]) and 
                                    zval >= sum(self.lengths_z[0:z]) and zval <= sum(self.lengths_z[0:z+1])):
                                    
                                    mesh.cells[zn*mesh.cells_y*mesh.cells_x + yn*mesh.cells_x + xn].setMaterial(cell.material)
                                    flag = 1
                                    break
                        
                            if flag == 1:
                                break

                        if flag == 1:
                            break

        return mesh


    # set teh boundaries for each cell
    def setCellBoundaries(self):

        # set x boundaries
        for z in range(self.cells_z):
            for y in range(self.cells_y):
                cell_num = z*self.cells_x*self.cells_y+y*self.cells_x
                self.cells[cell_num].surfaces[0].boundary = self.boundaries[0]
                cell_num = z*self.cells_x*self.cells_y+y*self.cells_x + self.cells_x-1
                self.cells[cell_num].surfaces[1].boundary = self.boundaries[1]


        # set y boundaries
        for z in range(self.cells_z):
            for x in range(self.cells_x):
                cell_num = z*self.cells_x*self.cells_y + x
                self.cells[cell_num].surfaces[2].boundary = self.boundaries[2]
                cell_num = z*self.cells_x*self.cells_y + (self.cells_y-1)*self.cells_x + x
                self.cells[cell_num].surfaces[3].boundary = self.boundaries[3]


        # set z boundaries
        for y in range(self.cells_y):
            for x in range(self.cells_x):
                cell_num = y*self.cells_x + x
                self.cells[cell_num].surfaces[4].boundary = self.boundaries[4]
                cell_num = (self.cells_z-1)*self.cells_x*self.cells_y + y*self.cells_x + x
                self.cells[cell_num].surfaces[5].boundary = self.boundaries[5]
