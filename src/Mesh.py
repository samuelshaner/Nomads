
from math import *
import numpy as np
from Cell import *
from Surface import *
import sys


class Mesh(object):

    def __init__(self, widths, heights):
        
        self.cells_x = len(widths)
        self.cells_y = len(heights)
        self.num_cells = self.cells_x * self.cells_y
        self.cells = np.empty(self.cells_x*self.cells_y, dtype=object)
        
        self.widths = widths
        self.heights = heights
        
        self.width = sum(widths)
        self.height = sum(heights)
        
        # create cell objects
        for y in range(self.cells_y):
            for x in range(self.cells_x):
                self.cells[y*self.cells_x+x] = Cell(widths[x], heights[y])
                
                
        # set neighbor cells and allocate for each cell
        for y in range(self.cells_y):
            for x in range(self.cells_x):
         
                cell = self.cells[y*self.cells_x+x]
         
                if x != 0:
                    cell.neighborCells[0] = self.cells[y*self.cells_x + x - 1]
        
                if y != self.cells_y - 1:
                    cell.neighborCells[1] = self.cells[(y+1)*self.cells_x + x]
                    
                if x != self.cells_x - 1:
                    cell.neighborCells[2] = self.cells[y*self.cells_x + x + 1]

                if y != 0:
                    cell.neighborCells[3] = self.cells[(y-1)*self.cells_x + x]


    def makeSurfaces(self):
        # set neighbor cells and allocate for each cell
        for y in range(self.cells_y):
            for x in range(self.cells_x):
         
                cell = self.cells[y*self.cells_x+x]
        
                # Surface 0
                if x == 0:
                    cell.surfaces[0] = Surface(cell.material.num_groups)
                else:
                    cell.surfaces[0] = self.cells[y*self.cells_x+x-1].surfaces[2]
                    
                # Surface 1
                cell.surfaces[1] = Surface(cell.material.num_groups)
                
                # Surface 2 
                cell.surfaces[2] = Surface(cell.material.num_groups)
                
                if y == 0:
                    cell.surfaces[3] = Surface(cell.material.num_groups)
                else:
                    cell.surfaces[3] = self.cells[(y-1)*self.cells_x+x].surfaces[1]
                    
    
    # break the geometry up into very small cells and solve diffusion problem
    def refineMesh(self, cell_size):
        
        # break up the old widths and heights into chunks .1 cm in size
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
                
                
        mesh = Mesh(widths, heights)
        
        for xn in range(mesh.cells_x):
            for yn in range(mesh.cells_y):
            
                xval = (sum(mesh.widths[0:xn]) + sum(mesh.widths[0:xn+1]))/2.0
                yval = (sum(mesh.heights[0:yn]) + sum(mesh.heights[0:yn+1]))/2.0
           
                flag = 0     
                for x in range(len(self.widths)):
                    for y in range(len(self.heights)):
                        
                        cell = self.cells[y*self.cells_x+x]
                        
                        if xval >= sum(self.widths[0:x]) and xval <= sum(self.widths[0:x+1]) and yval >= sum(self.heights[0:y]) and yval <= sum(self.heights[0:y+1]):
                            mesh.cells[yn*mesh.cells_x + xn].setMaterial(cell.material, cell.order)
                            flag = 1
                            break
                        
                    if flag == 1:
                        break
                
        return mesh
                
        
                    
                    
                    
                    
                    