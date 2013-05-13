
from math import *
import numpy as np

class Material(object):
    
    def __init__(self, num_groups, id):
        self.D           = np.zeros(num_groups)
        self.num_groups  = num_groups
        self.sigma_a     = np.zeros(num_groups)
        self.sigma_t     = np.zeros(num_groups)
        self.nu_sigma_f  = np.zeros(num_groups)
        self.sigma_s     = np.zeros([num_groups, num_groups])
        self.chi         = np.zeros(num_groups)
        self.id          = id

    def setSigmaS(self, xs):
        
        self.sigma_s = xs
    
    def setSigmaT(self, xs):
        
        self.sigma_t = xs
        
    def setSigmaA(self, xs):
        
        self.sigma_a = xs
        
    def setNuSigmaF(self, xs):
        
        self.nu_sigma_f = xs
        
    def setChi(self, chi):
        
        self.chi = chi
        
    def setD(self, D):
        
        self.D = D
                
        
        
        
        
        
        
        