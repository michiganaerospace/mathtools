'''An interpolation object for fitting data with Splines, Legendre Polynomials, 
and Fourier Series.'''
import numpy as np    
from mathtools.legendre import *
from mathtools.utils import Struct, validate_type
import pdb


# Fit class. A general purpose regularized interpolation machine.
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

        self.basis_generators = {}
        self.basis_generators['legendre'] = create_legendre_basis
        self.coefs = None

        # Assign variables to the object.
        self.x = x
        self.basis_type = basis_type
        self.nb_bases = nb_bases
        self.reg_coefs = reg_coefs

        # Error checking.
        self._validate()

        # We're clean. Let's generate a basis!
        self.basis = self.basis_generators[basis_type](x, nb_bases, reg_coefs)


    def config(self, x=None, nb_bases=None, basis_type=None, reg_coefs=None):
        ''' Reconfigure the fit object with new bases, regularization
            coefficients, etc.
        '''
        if (x is not None):
            self.x = x
        if nb_bases is not None:
            self.nb_bases = nb_bases
        if (basis_type is not None):
            self.basis_type = basis_type
        if (reg_coefs is not None):
            self.reg_coefs = reg_coefs

        # Now recompute the basis object with these changes.
        self.basis = self.basis_generators[self.basis_type]\
                (self.x, self.nb_bases, self.reg_coefs)

        # Ensure that fit coefficients are set to zero.
        self.coefs = None


    def fit(self, y):
        ''' Find the regularized least squares fit to the data y.'''

        # Error checking!
        if len(y) != len(self.x):
            raise ValueError('Data is not the same length as domain!')

        # Fit the data!
        fit = least_squares(self.basis, y)
        self.coefs = fit.coefs
        self.current_fit = fit

        return fit


    def resample(self, x):
        '''Resample the current fit to a different domain.'''

        # We must have fit data before we can resample.
        if (self.coefs is None):
            raise ValueError('No data has been fit. Cannot resample.')
        
        # Generate a new basis.
        resampled_basis = self.basis_generators[self.basis_type]\
                (x, self.nb_bases, x_ref=self.x)

        # Project coefficients onto the resampled basis.
        fit = least_squares(resampled_basis, coefs=self.coefs)

        return fit
    

    def _validate(self):
        # Validate the parameters of the fit object. 

        # Validate the types.
        validate_type(self.x, [list, np.ndarray], 'x')
        validate_type(self.nb_bases, [int], 'nb_bases')

        # Check values.
        if self.nb_bases == 0:
            raise ValueError('You must specify a nonzero number of basis\
                              vectors')

        if self.basis_type not in self.available_bases:
            raise ValueError('{} is not a valid basis type.'\
                    .format(self.basis_type))


def least_squares(basis, y=None, coefs=None):
    '''Find the least square fit to one-dimensional data using the specified 
       basis. OR, if coefficients are specified rather than data, y, the
       coefficients are used to generate the fit; in this case, the number of
       coefficients must equal basis.nb_bases.
    INPUT
        basis - basis object
            A mathtools basis object (see mathtools.legendre, e.g.)
        y - array_like
            One dimensional array of data to fit.
        coefs - array_like
            Coefficients used to project basis generate fit.
    OUTPUT
       fit - Struct
            A Struct containing the following fields:
                - x:    The domain on which fit is defined.
                - y:    The best fit to the data.
                - dy:   The derivative of the best fit.
                - d2y:  The second derivative of the best fit.
                - coefs The coefficients used for the fit.
    '''

    # Do we have data, or must we do the fit?
    if (y is not None):
        # Augment the data (for regularization).
        augmented_y = basis.augment(y)

        # Perform the least squares fit using the precomputed pseudoinverse!
        coefs = basis.inverse.dot(augmented_y)

    if (coefs is None):
        raise ValueError('Cannot project onto basis! No data or coefficients!')

    # Generate a fit Struct object to hold the results.
    fit                         = Struct()
    fit.coefs                   = coefs
    fit.x                       = basis.x
    fit.y                       = np.zeros_like(basis.x)
    fit.dy                      = np.zeros_like(basis.x)
    fit.d2y                     = np.zeros_like(basis.x)
    fit.y[basis.valid_idx]      = basis.B.dot(coefs)
    fit.dy[basis.valid_idx]     = basis.dB.dot(coefs)
    fit.d2y[basis.valid_idx]    = basis.d2B.dot(coefs)

    return fit

