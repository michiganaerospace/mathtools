from mathtools.fit import *
from nose.tools import assert_equals, assert_almost_equals, assert_raises
from nose import with_setup
from numpy.testing import assert_array_almost_equal_nulp, assert_array_equal
import numpy as np


# SETUP -----------------------------------------------------

def setup():
    global t, y
    t = np.linspace(0,5*np.pi,200)
    y = np.sin(2*np.pi/5*t) + 0.2 * np.random.randn(len(t))


def teardown():
    pass


# BEGIN TESTS ------------------------------------------------------

def create_struct_test():
    s = Struct()
    s.name = 'Matthew'
    s.id = 31415
    assert_equals(s.name, 'Matthew')
    assert_equals(s.id, 31415)


@with_setup(setup, teardown)
def create_unknown_basis_test():
    assert_raises(ValueError, create_basis, t, 5, 'Matthew')


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
    # Size of brick is 15 x (3*len(t) + 15)
    f = Fit(t, y, 15)
    assert_equals(f.inverse.shape, (15,3*len(t) + 15))


@with_setup(setup, teardown)
def augmented_data_test():
    lt = len(t)
    nb_orders = 15
    f = Fit(t, y, nb_orders)
    assert_equals(f.y_aug.shape, (3*lt+nb_orders,))


@with_setup(setup, teardown)
def precompute_bases_test():
    f = Fit(x=t, nb_orders=15)
    assert_equals(f.inverse.shape, (15,615))


@with_setup(setup, teardown)
def bases_present_test():
    f = Fit(x=t, nb_orders=15)
    assert_equals(f.B.shape, (200,15))


@with_setup(setup, teardown)
def identity_present_test():
    f = Fit(x=t, nb_orders=15)
    assert_equals(f.I.shape, (15,15))


@with_setup(setup, teardown)
def derivative_present_test():
    f = Fit(x=t, nb_orders=15)
    assert_equals(f.dB.shape, (200,15))


@with_setup(setup, teardown)
def second_derivative_present_test():
    f = Fit(x=t, nb_orders=15)
    assert_equals(f.d2B.shape, (200,15))


@with_setup(setup, teardown)
def raises_error_on_domain_mismatch_test():
    bad_t = np.linspace(0,5*np.pi,100)
    assert_raises(ValueError, Fit, bad_t, y, 15)


@with_setup(setup, teardown)
def bases_struct_test():
    f = Fit(t, y, nb_orders=15, reg_coefs=[0,1e-3,1e-3])
    assert_array_equal(f.bases.B, f.B)
    assert_array_equal(f.bases.dB, f.dB)
    assert_array_equal(f.bases.d2B, f.d2B)
    assert_array_equal(f.bases.B_, f.B_)
    assert_equals(str(type(f.bases)), "<class 'mathtools.utils.Struct'>") 


@with_setup(setup, teardown)
def results_struct_test():
    f = Fit(t, y, nb_orders=15, reg_coefs=[0,1e-3,1e-3])
    assert_array_equal(f.results.x, f.x)
    assert_equals(len(f.results.coefs), 15)
    assert_equals(f.results.y.shape, f.x.shape)
    assert_equals(f.results.dy.shape, f.x.shape)
    assert_equals(f.results.d2y.shape, f.x.shape)

#TODO -- Implement resampling.

if __name__ == "__main__":
    from mathtools.vanity import *
    setup_plotting() 

    t = np.linspace(0,5*np.pi, 300)
    y = np.sin(2*np.pi/5*t) + 0.2 * np.random.randn(len(t))
    f = Fit(t,y, 15, reg_coefs=[0,1e-3,1e-4])

    figure(100)
    plot(t, y, color=pomegranate, linewidth=2, alpha=0.6)
    hold(True)
    plot(f.results.x, f.results.y, color=belize_hole, linewidth=2)
    grid(True)
    
