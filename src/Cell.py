
from math import *
import numpy as np
from Surface import *

class Cell(object):
    
    def __init__(self, width, length, height):
        self.surfaces = np.empty(6, dtype=object)
        self.neighborCells = np.empty(6, dtype=object)    
        self.width = width
        self.length = length
        self.height = height
        self.volume = width*length*height        
        
    def setMaterial(self, material):
        self.material = material
        self.flux     = np.ones(material.num_groups)
        self.current  = np.zeros(material.num_groups)
        self.old_flux = np.ones(material.num_groups)
