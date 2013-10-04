
from math import *
import numpy as np
from Surface import *

class Cell(object):
    
    def __init__(self, length_x, length_y, length_z):
        self.surfaces = np.empty(6, dtype=object)
        self.neighborCells = np.empty(6, dtype=object)    
        self.length_x = length_x
        self.length_y = length_y
        self.length_z = length_z
        self.volume   = length_x*length_y*length_z        
        
    def setMaterial(self, material):
        self.material = material
        self.flux     = np.ones(material.num_groups)
        self.current  = np.zeros(material.num_groups)
        self.old_flux = np.ones(material.num_groups)
