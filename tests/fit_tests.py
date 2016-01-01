from mathtools.fit import *
from nose.tools import assert_equals, assert_almost_equals, assert_raises
from nose import with_setup
from numpy.testing import assert_array_almost_equal_nulp, assert_array_equal
import numpy as np


# BEGIN TESTS ------------------------------------------------------

# Can you create a fit object?
def init_test():
    f = Fit()
    assert_equals(str(type(f)), "<class 'mathtools.fit.Fit'>") 


def legendre_by_default_test():
    f = Fit()
    assert_equals(f.basis_type, 'legendre') 


def default_regularization_coefficients_test():
    f = Fit()
    assert_array_almost_equal_nulp(f.reg_coefs,[0.0, 0.0, 0.0]) 


def legendre_basis_test():
    t = np.linspace(0,2*np.pi)
    y = np.sin(2*np.pi/5*t) + 0.2 * np.random.randn(len(t))
    f = Fit(t, y, 15)
    assert_equals(f.B.shape, (len(t), 15))
    assert_equals(f.dB.shape, (len(t), 15))
    assert_equals(f.d2B.shape, (len(t), 15))


def legendre_basis_test():
    t = np.linspace(0,2*np.pi)
    y = np.sin(2*np.pi/5*t) + 0.2 * np.random.randn(len(t))
    f = Fit(t, y, 15)
    assert_equals(f.B.shape, (len(t), 15))
    assert_equals(f.dB.shape, (len(t), 15))
    assert_equals(f.d2B.shape, (len(t), 15))

if __name__ == '__main__':
    
    t = np.linspace(0,2*np.pi)
    y = np.sin(2*np.pi/5*t) + 0.2 * np.random.randn(len(t))
    f = Fit(t, y, nb_orders=15)
