
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
                    
                    
                    
                    
                    