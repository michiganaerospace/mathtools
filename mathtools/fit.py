'''An interpolation object for fitting data with Splines, Legendre Polynomials, 
and Fourier Series.'''
import numpy as np    
from mathtools.bases import *
import pdb

# Fit class. A general purpose interpolation machine.
class Fit(object):

    def __init__(self, x=None, y=None, nb_orders=0, basis_type='legendre',  \
                 reg_coefs=[0.0, 0.0, 0.0], existing_basis=None,            \
                 filename=None):

        # Assign object properties.
        self.basis_type = basis_type
        self.x = x
        self.y = y
        self.nb_orders = nb_orders
        self.reg_coefs = np.array(reg_coefs)

        if (self.x is not None) and (self.nb_orders>0):
            self._create_basis()
        
    def _create_basis(self):
        # Call appropriate basis generation.
        if self.basis_type == 'legendre':
           self._create_legendre_basis() 
        elif self.basis_type == 'fourier':
            pass
        elif self.basis_type == 'cubic-spline':
            pass

    def _create_legendre_basis(self):
        # Create the legendre bases.
        self.B   = legendre_basis(self.x, self.nb_orders)
        self.dB  = legendre_basis(self.x, self.nb_orders)
        self.d2B = legendre_basis(self.x, self.nb_orders)


       
