'''An interpolation object for fitting data with Splines, Legendre Polynomials, 
and Fourier Series.'''
import numpy as np    
from mathtools.bases import *
import pdb

class Details(object):
    pass


# Fit class. A general purpose interpolation machine.
class Fit(object):


    def __init__(self, x=None, y=None, nb_orders=0, basis_type='legendre',  \
                 reg_coefs=[0.0, 0.0, 0.0])

        # Assign object properties.
        self.basis_type = basis_type
        self.x          = x
        self.y          = y
        self.nb_orders  = nb_orders
        self.reg_coefs  = np.array(reg_coefs)
        self.inverse    = None 
        self.fits       = Details()

        # Create the fitting basis.
        if (self.x is not None) and (self.nb_orders>0):
            self._create_basis()

        # IF we have enough data, fit the data.
        if (self.y is not None) and (self.inverse is not None):
            self._augment_y()
            self._compute_fit_coefficients()
            self._fit() 

    def __call__(self):
        print 'Called!'
            
    def _create_basis(self):
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
        # Create the legendre bases.
        reg_coefs   = self.reg_coefs
        self.B      = legendre_basis(self.x, self.nb_orders)
        self.I      = reg_coefs[0] * np.eye(self.nb_orders)
        self.dB     = reg_coefs[1] * legendre_basis(self.x, self.nb_orders)
        self.d2B    = reg_coefs[2] * legendre_basis(self.x, self.nb_orders)

        # Create the 'brick' by stacking these bases on top of one another.
        self.B_     = np.r_[self.B, self.I, self.dB, self.d2B] 

        # Compute the inverse operator for the least squares problem.
        self._compute_inverse()

    def _compute_inverse(self):
        # Find the inverse associated with the 'brick'. Keep it around for
        # computational efficiency. Using the pinv function, which in turn
        # uses numpy's SVD to compute the pseudoinverse of the brick B_.
        self.inverse = np.linalg.pinv(self.B_)

    def _augment_y(self):
        # Augment the data vectors with appropriate zeros in order to allow for
        # regularization on coefficient magnitudes and on the magnitude of the
        # derivatives.
        self.y_aug = np.r_[self.y, np.zeros(2*len(self.x)+self.nb_orders)]

    def _compute_fit_coefficients(self):
        # Compute the regularized fit coefficients.
        self.fit_coefs = self.inverse.dot(self.y_aug)
        
    def _fit(self):
        self.fits.x      = self.x
        self.fits.coefs  = self.fit_coefs
        self.fits.y      = self.B.dot(self.fit_coefs)
        self.fits.dy     = self.dB.dot(self.fit_coefs)
        self.fits.d2y    = self.d2B.dot(self.fit_coefs)
   
