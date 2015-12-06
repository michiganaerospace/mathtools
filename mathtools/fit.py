'''An interpolation object for fitting data with Splines, Legendre Polynomials, 
and Fourier Series.'''
import numpy as np    


# Fit class. A general purpose interpolation machine.
class Fit(object):

    def __init__(self, x=None, y=None, nb_orders=0, basis_type='legendre',  \
                 reg_y=0, reg_dy=0, reg_d2y=0, existing_basis=None,         \
                 filename=None):

        # Assign object properties.
        self.basis_type = basis_type

    def _create_basis(self):
        # Create the basis vectors on the domain.
        pass
