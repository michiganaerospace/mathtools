'''An interpolation object for fitting data with Splines, Legendre Polynomials, 
and Fourier Series.'''
import numpy as np    
from mathtools.legendre import *
from mathtools.utils import Struct
import pdb


def pseudoinverse(M):
    '''Find the pseudoinverse of the matrix M using singular value
       decomposition.
    INPUT
        M - array_like
            An (m x n) matrix whose pseudoinverse we want to find.
    OUTPUT
    '''
    U,s,V_T = np.linalg.svd(M)


def create_basis(x, nb_bases, basis_type='legendre', x_reference=None):
    '''Create a basis structure.
    INPUTS
        x - array_like
            Domain over which we wish to build the basis matrix.
        nb_bases - int
            The number of basis vectors to use.
        basis_type - str
            Type of basis functions to uses. May be 'legendre', 'fourier', or
            'cubic-spline'.
    OUTPUT
        basis - Struct object
            A Struct object containing the following properties:
                - B
                - dB
                - d2B
                - I
                - 
    '''
    # Error checking.
    if basis_type not in ['legendre', 'fourier', 'cubic-spline']:
        raise ValueError('Unknown basis type requested.')

    # Create hold structure.
    basis = Struct()

    # Build basis.

# Fit class. A general purpose interpolation machine.
class Fit(object):

    def __init__(self, x=None, y=None, nb_orders=0, basis_type='legendre',\
            reg_coefs=[0.0, 0.0, 0.0]):

        # Assign object properties.
        self.basis_type = basis_type
        self.x          = x
        self.y          = y
        self.nb_orders  = nb_orders
        self.reg_coefs  = np.array(reg_coefs)
        self.inverse    = None
        self._init()

    def _init(self):
        # Initialize everything we can initialize, given what has been
        # provided.
        self.inverse    = None 
        self._results   = Struct()

        # Create the basis matrix.
        if self._can_compute_basis():
            self._create_basis()

        # IF we have enough data, fit the data.
        if self._can_fit():
            self._fit()
        
    def config(self, x=None, y=None, nb_orders=None, basis_type=None,\
            reg_coefs=None):
        # Configure the current Fit object to have new abscissa, order,
        # basis-type, etc.
        if (x is not None):
            self.x = x
        if (y is not None):
            self.y = y
        if (nb_orders is not None):
            self.nb_orders = nb_orders
        if (basis_type is not None):
            self.basis_type = basis_type
        if (reg_coefs is not None):
            self.reg_coefs = reg_coefs

        # Recreate bases, etc., as necessary.
        self._init()

    def _fit(self):
        # Function that uses precomputed bases to find regularized least 
        # squares fit to the currently observed data, y.
        if self._can_fit():
            self._augment_y()
            self._compute_fit_coefficients()
            self._compute_fit() 

    def fit(self, x=None, y=None):
        # Fit the data y. If new abscissa x is provided, recompute the bases.
        # Return the data results object to the user.

        # If x is provided, recompute the bases.
        if (x is not None):
            self.x = x
            self._create_basis()

        # If a y vector is provided, fit the data!
        if (y is not None):
            self.y = y

        # Do the fitting to the data.
        self._fit()
        
        return self._results

    def _can_fit(self):
        return (self.x is not None) and (self.y is not None) and \
               (self.nb_orders > 1) and \
               (self.basis_type in ['legendre', 'cubic-spline', 'fourier'])

    def _can_compute_basis(self):
        return (self.x is not None) and (self.nb_orders > 1)

    def _build_basis(self, x, basis_type=None, x_orig=None):
        # Build a basis of specified type given the data points x. If shift and
        # scale numbers are provided, data will be shifted and scaled -- x_ =
        # scale * (x_ - shift) before the basis is created. 
        pass
        

    def _create_basis(self):
        # Compute the basis vector matrices.

        # Call appropriate basis generation.
        if self.basis_type == 'legendre':
           self._create_legendre_basis() 
        elif self.basis_type == 'fourier':
            # TODO
            pass
        elif self.basis_type == 'cubic-spline':
            # TODO
            pass

    def _create_legendre_basis(self):
        # Create the legendre bases for y and its derivatives.

        # Build bases and the 'brick'.
        reg_coefs   = self.reg_coefs
        self.B      = legendre_basis(self.x, self.nb_orders)
        self.I      = reg_coefs[0] * np.eye(self.nb_orders)
        self.dB     = reg_coefs[1] * d_legendre_basis(self.x, self.nb_orders)
        self.d2B    = reg_coefs[2] * d2_legendre_basis(self.x, self.nb_orders)

        # Create the 'brick' by stacking these bases on top of one another.
        self.B_     = np.r_[self.B, self.I, self.dB, self.d2B] 

        # Compute the inverse operator for the least squares problem.
        self._compute_inverse()

    def _compute_inverse(self):
        # Find the inverse associated with the 'brick'. Keep it around for
        # computational efficiency. Using the pinv function, which in turn
        # uses numpy's SVD to compute the pseudoinverse of the brick B_.
        U,s,V_transpose = np.linalg.svd(self.B_)
        self.inverse = np.linalg.pinv(self.B_)

    def _augment_y(self):
        # Augment the data vectors with appropriate zeros in order to allow for
        # regularization on coefficient magnitudes and on the magnitude of the
        # derivatives.
        self.y_aug = np.r_[self.y, np.zeros(2*len(self.x)+self.nb_orders)]

    def _compute_fit_coefficients(self):
        # Compute the regularized fit coefficients.
        if (len(self.x) != len(self.y)):
            raise ValueError('Vectors x and y must have the same length!')
        self.fit_coefs = self.inverse.dot(self.y_aug)
        
    def _compute_fit(self):
        self._results.x      = self.x
        self._results.coefs  = self.fit_coefs
        self._results.y      = self.B.dot(self.fit_coefs)
        self._results.dy     = self.dB.dot(self.fit_coefs)
        self._results.d2y    = self.d2B.dot(self.fit_coefs)

    def resample(self, x=None):
        # Resample the current fit to a different domain.
        self.x_resampled = x
        #TODO -- Resample in some smart way here. Probably need refactor.
   
    @property
    def results(self):
        return self._results

    @property
    def bases(self):
        # For convenience, package current bases into a single
        # structure and send it out to the user.
        b       = Struct()
        b.B     = self.B
        b.dB    = self.dB
        b.d2B   = self.d2B
        b.B_    = self.B_
        return b
