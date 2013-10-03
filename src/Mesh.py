
from math import *
import numpy as np
from Cell import *
from Surface import *
import sys


class Mesh(object):

    def __init__(self, widths, lengths, heights=[]):
        
        self.cells_x = len(widths)
        self.cells_y = len(lengths)
        
        self.widths = widths
        self.lengths = lengths
        
        self.width = sum(widths)
        self.length = sum(lengths)

        self.boundaries = np.empty(6, dtype=object)
        for i in range(6):
            self.boundaries[i] = Surface(0)

        if heights:
            self.dimensions = 3
            self.heights = heights
            self.height = sum(heights)
            self.cells_z = len(heights)
        else:
            self.dimensions = 2
            self.heights = [1.0]
            self.height = 1.0
            self.cells_z = 1

        self.num_cells = self.cells_x * self.cells_y * self.cells_z
        self.cells = np.empty(self.num_cells, dtype=object)

        # create cell objects
        for z in range(self.cells_z):
            for y in range(self.cells_y):
                for x in range(self.cells_x):
                    self.cells[z*(self.cells_x*self.cells_y) + y*self.cells_x + x] = Cell(self.widths[x], self.lengths[y], self.heights[z])
                                
        # set neighbor cells and allocate for each cell
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
                    

    def setBoundaries(self, left, right, back, front, bottom='REFLECTIVE', top='REFLECTIVE'):

        self.boundaries[0].boundary = left
        self.boundaries[1].boundary = right
        self.boundaries[2].boundary = back
        self.boundaries[3].boundary = front
        self.boundaries[4].boundary = bottom
        self.boundaries[5].boundary = top


    def makeSurfaces(self):
        # set neighbor cells and allocate for each cell
        for z in range(self.cells_z):
            for y in range(self.cells_y):
                for x in range(self.cells_x):
         
                    cell_num = z*self.cells_x*self.cells_y + y*self.cells_x + x
                    cell = self.cells[cell_num]
        
                    # Surface 0
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

                    # Surface 1
                    cell.surfaces[1] = Surface(cell.material.num_groups)
                
                    # Surface 3
                    cell.surfaces[3] = Surface(cell.material.num_groups)

                    # Surface 5
                    cell.surfaces[5] = Surface(cell.material.num_groups)                
                    
    
    # break the geometry up into very small cells and solve diffusion problem
    def refineMesh(self, cell_size, dimensions=[0,1]):
        
        if 0 in dimensions:
            # break up the old widths into chunks
            widths = []
            for i in self.widths:
                num_new_widths = int(floor(i/cell_size))
                new_widths = []

                for j in range(num_new_widths):
                    new_widths.append(cell_size)
                
                if abs(sum(new_widths) - i) >= 1.e-3:
                    new_widths.append(i - sum(new_widths))
                
                for k in new_widths:
                    widths.append(k) 
        else:
            widths = self.widths

        if 1 in dimensions:
            # break up the old lengths into chunks
            lengths = []
            for i in self.lengths:
                num_new_lengths = int(floor(i/cell_size))
                new_lengths = []
            
                for j in range(num_new_lengths):
                    new_lengths.append(cell_size)
                
                if abs(sum(new_lengths) - i) >= 1.e-3:
                    new_lengths.append(i - sum(new_lengths))
                
                for k in new_lengths:
                    lengths.append(k)

        else:
            lengths = self.lengths

        if 2 in dimensions:
            # break up the old heights into chunks
            heights = []
            for i in self.heights:
                num_new_heights = int(floor(i/cell_size))
                new_heights = []
            
                for j in range(num_new_heights):
                    new_heights.append(cell_size)
                
                if abs(sum(new_heights) - i) >= 1.e-3:
                    new_heights.append(i - sum(new_heights))
                
                for k in new_heights:
                    heights.append(k)
        else:
            heights = self.heights

        if self.dimensions == 2:
            mesh = Mesh(widths, lengths)
        else:
            mesh = Mesh(widths, lengths, heights)

        # set the cell materials
        for zn in range(mesh.cells_z):
            for yn in range(mesh.cells_y):
                for xn in range(mesh.cells_x):
            
                    xval = (sum(mesh.widths[0:xn]) + sum(mesh.widths[0:xn+1]))/2.0
                    yval = (sum(mesh.lengths[0:yn]) + sum(mesh.lengths[0:yn+1]))/2.0
                    zval = (sum(mesh.heights[0:zn]) + sum(mesh.heights[0:zn+1]))/2.0
           
                    flag = 0     
                    for x in range(len(self.widths)):
                        for y in range(len(self.lengths)):
                            for z in range(len(self.heights)):
                                
                                cell_num = z*self.cells_x*self.cells_y + y*self.cells_x + x
                                cell = self.cells[cell_num]
                        
                                if (xval >= sum(self.widths[0:x]) and xval <= sum(self.widths[0:x+1]) and yval >= sum(self.lengths[0:y]) and yval <= sum(self.lengths[0:y+1]) 
                                    and zval >= sum(self.heights[0:z]) and zval <= sum(self.heights[0:z+1])):
                                    mesh.cells[zn*mesh.cells_y*mesh.cells_x + yn*mesh.cells_x + xn].setMaterial(cell.material)
                                    flag = 1
                                    break
                        
                            if flag == 1:
                                break

                        if flag == 1:
                            break

        return mesh


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



