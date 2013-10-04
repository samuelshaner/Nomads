
from math import *
import numpy as np

class Material(object):
    
    def __init__(self, num_groups, id):
        
        # initialize the cross section arrays
        self.D           = np.zeros(num_groups)
        self.sigma_a     = np.zeros(num_groups)
        self.sigma_r     = np.zeros(num_groups)
        self.sigma_t     = np.zeros(num_groups)
        self.nu_sigma_f  = np.zeros(num_groups)
        self.sigma_s     = np.zeros([num_groups, num_groups])
        self.chi         = np.zeros(num_groups)

        # set other material properties
        self.id          = id
        self.buckling    = 0.0
        self.num_groups  = num_groups


    # set the scattering cross sections
    def setSigmaS(self, xs):
        
        self.sigma_s = xs

        for e in range(self.num_groups):
            for g in range(self.num_groups):
                if g != e:
                    self.sigma_r[e] += self.sigma_s[e,g]

    # set the total cross sections
    def setSigmaT(self, xs):
        
        self.sigma_t = xs
        
    # set the absorption cross sections
    def setSigmaA(self, xs):
        
        self.sigma_a = xs

        for e in range(self.num_groups):
            self.sigma_r[e] += self.sigma_a[e]

    # set the nu sigma fission cross sections
    def setNuSigmaF(self, xs):
        
        self.nu_sigma_f = xs
        
    # set the fission spectrum
    def setChi(self, chi):
        
        self.chi = chi

    # set the buckling coefficient and increment the absorption and removal xs
    def setBuckling(self, buckling):
        
        self.buckling = buckling

        if self.D[0] > 0.0:
            for e in range(self.num_groups):
                self.sigma_r[e] += self.D[e] * self.buckling
                self.sigma_a[e] += self.D[e] * self.buckling
        
    # set the diffusion coefficients and increment the absorption and removal xs
    def setD(self, D):
        
        self.D = D

        if self.buckling > 0.0:
            for e in range(self.num_groups):
                self.sigma_r[e] += self.D[e] * self.buckling
                self.sigma_a[e] += self.D[e] * self.buckling        

        
