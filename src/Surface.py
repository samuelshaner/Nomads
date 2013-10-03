
from math import *
import numpy as np

class Surface(object):
    
    def __init__(self, num_groups):
        self.current  = np.zeros(num_groups)
        self.DDif     = np.zeros(num_groups)
        self.DHat     = np.zeros(num_groups)
        self.boundary = 'REFLECTIVE'
        
