'''An interpolation object for fitting data with Splines, Legendre Polynomials, 
and Fourier Series.'''
import numpy as np    
from mathtools.legendre import *
from mathtools.utils import Struct
import pdb


# Fit class. A general purpose interpolation machine.
class Fit(object):

    def __init__(self, x, nb_bases=0, basis_type='legendre', reg_coefs=[0,0,0]):
        '''Initialize the fit object. 
        INPUTS
            x - array_like
                The domain on which we would like to fit things.
            nb_bases - int
                Number of basis vectors to use.
            basis_type - string [default: 'legendre']
                The type of basis to use. Must be of type 'legendre'.
            reg_coefs - array_like
                A length 3 array of of regularization coefficients that
                determine how much to penalize the solution magnitude, the
                magnitude of the derivative and magnitude of the second
                derivative, respectively.
        '''

        # Currently available bases.
        self.available_bases = ['legendre']

        self.basis_machines = {}
        self.basis_machines['legendre'] = create_legendre_basis

        # Error checking.
        if nb_bases == 0:
            raise ValueError('You must specify a nonzero number of basis\
                              vectors')

        if basis_type not in self.available_bases:
            raise ValueError('{} is not a valid basis type.'.format(basis_type))

        # We're clean. Assign variables to object.
        self.basis_type = basis_type
        self.basis = self.basis_machines[basis_type](x, nb_bases, reg_coefs)
        



        
