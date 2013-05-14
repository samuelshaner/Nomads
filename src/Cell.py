
from math import *
import numpy as np
from Surface import *

class Cell(object):
    
    def __init__(self, width, height):
        self.surfaces = np.empty(4, dtype=object)
        self.neighborCells = np.empty(4, dtype=object)    
        self.width = width
        self.height = height
        self.volume = width*height        
        
    def setMaterial(self, material, order):
        self.material = material
        self.flux     = np.ones(material.num_groups)
        self.current  = np.zeros(material.num_groups)
        self.old_flux = np.ones(material.num_groups)
        self.coeffs   = np.zeros(order*material.num_groups)
    