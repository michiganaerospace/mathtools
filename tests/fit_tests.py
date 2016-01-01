from mathtools.fit import *
from nose.tools import assert_equals, assert_almost_equals, assert_raises
from nose import with_setup
from numpy.testing import assert_array_almost_equal_nulp, assert_array_equal
import numpy as np


# SETUP -----------------------------------------------------

def setup():
    global t, y
    t = np.linspace(0,2*np.pi)
    y = np.sin(2*np.pi/5*t) + 0.2 * np.random.randn(len(t))


def teardown():
    pass


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
    assert_array_almost_equal_nulp(f.reg_coefs, [0.0, 0.0, 0.0]) 


@with_setup(setup, teardown)
def legendre_basis_test():
    f = Fit(t, y, 15)
    assert_array_almost_equal_nulp(f.x, t)
    assert_array_almost_equal_nulp(f.y, y)
    assert_equals(f.B.shape, (len(t), 15))
    assert_equals(f.dB.shape, (len(t), 15))
    assert_equals(f.d2B.shape, (len(t), 15))


@with_setup(setup, teardown)
def brick_creation_test():
   f = Fit(t, y, 15) 
   assert_equals(f.B_.shape, (3*len(t)+15, 15))


@with_setup(setup, teardown)
def inverse_creation_test():
    # size of brick is 15 x (3*len(t) + 15)
    f = Fit(t, y, 15)
    assert_equals(f.inverse.shape, (15,150+15))


@with_setup(setup, teardown)
def augmented_data_test():
    lt = len(t)
    nb_orders = 15
    f = Fit(t, y, nb_orders)
    assert_equals(f.y_aug.shape, (3*lt+nb_orders,))


if __name__ == "__main__":
    from mathtools.vanity import *
    setup_plotting() 

    t = np.linspace(0,5*np.pi, 300)
    y = np.sin(2*np.pi/5*t) + 0.2 * np.random.randn(len(t))
    f = Fit(t,y, 15, reg_coefs=[0,0.1,0.2])

    figure(100)
    plot(t, y, color=pomegranate, linewidth=2, alpha=0.6)
    hold(True)
    plot(f.fit.x, f.fit.y, color=belize_hole, linewidth=2)
    grid(True)
    
