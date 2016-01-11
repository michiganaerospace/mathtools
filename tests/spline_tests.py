
'''Test basis generation algorithms for Fourier bases.'''
from mathtools.fit import *
import mathtools.refs.FS_Tools as fs
from mathtools.utils import *
from mathtools.splines import *
from nose.tools import assert_equals, assert_almost_equals, assert_raises
from numpy.testing import assert_array_almost_equal_nulp, assert_array_equal
from nose import with_setup
import numpy as np


# BEGIN TESTS ------------------------------------------------------

# Create a spline basis? 
def create_spline_basis_on_unit_interval_test():
    x = np.linspace(0,5,100)
    B = cubic_spline_basis_unit_interval(x)
    assert_equals(B.shape, (100,4))


def unit_spline_matches_ref_test():
    x = np.linspace(0,5,100)
    B = cubic_spline_basis_unit_interval(x)
    B_fs = fs.cubic_spline_interval(x)
    assert_array_almost_equal_nulp(B, B_fs)


def create_spline_d_basis_on_unit_interval_test():
    x = np.linspace(0,5,100)
    dB = d_cubic_spline_basis_unit_interval(x)
    assert_equals(dB.shape, (100,4))


def d_unit_spline_matches_ref_test():
    x = np.linspace(0,5,100)
    dB = d_cubic_spline_basis_unit_interval(x)
    dB_fs = fs.dcubic_spline_interval(x)
    assert_array_almost_equal_nulp(dB, dB_fs)


def create_spline_d2_basis_on_unit_interval_test():
    x = np.linspace(0,5,100)
    d2B = d2_cubic_spline_basis_unit_interval(x)
    assert_equals(d2B.shape, (100,4))


def d2_unit_spline_matches_ref_test():
    x = np.linspace(0,5,100)
    d2B = d2_cubic_spline_basis_unit_interval(x)
    d2B_fs = fs.d2cubic_spline_interval(x)
    assert_array_almost_equal_nulp(d2B, d2B_fs)


def knot_generation_works_test():
    x = np.linspace(0,5,100)
    knots = uniform_knots(x, 20)
    knots_fs = fs.get_uniform_knots(x, 20)
    assert_array_almost_equal_nulp(knots, knots_fs)


def generate_full_cubic_spline_basis_test():
    x = np.linspace(0,5,100)
    knots = uniform_knots(x, 20)
    B = cubic_spline_basis(x, 20)
    B_fs = fs.get_cubic_spline_basis(x, knots)
    assert_array_almost_equal_nulp(B, B_fs)
    

def generate_full_d_cubic_spline_basis_test():
    x = np.linspace(0,5,100)
    knots = uniform_knots(x, 20)
    dB = d_cubic_spline_basis(x, 20)
    dB_fs = fs.get_dcubic_spline_basis(x, knots)
    assert_array_almost_equal_nulp(dB, dB_fs)


def generate_full_d2_cubic_spline_basis_test():
    x = np.linspace(0,5,100)
    knots = uniform_knots(x, 20)
    d2B = d2_cubic_spline_basis(x, 20)
    d2B_fs = fs.get_d2cubic_spline_basis(x, knots)
    assert_array_almost_equal_nulp(d2B, d2B_fs)


def create_spline_test_01():
    x = np.linspace(0,5,100)
    nb_bases = 25
    basis = create_spline_basis(x, nb_bases)
    assert_equals(basis.nb_bases, 25)


def create_spline_shape_test():
    x = np.linspace(0,5,100)
    nb_bases = 25
    basis = create_spline_basis(x, nb_bases)
    assert_equals(basis.B.shape, (100, 25))
    assert_equals(basis.dB.shape, (100, 25))
    assert_equals(basis.d2B.shape, (100, 25))
    assert_equals(basis.B_.shape, (325, 25))


def create_spline_shape_test():
    x = np.linspace(0,5,100)
    nb_bases = 25
    basis = create_spline_basis(x, nb_bases, reg_coefs=[1e-3, 1e-3, 1e-3])
    assert_equals(basis.inverse.shape, (25, 325))


def create_spline_condition_number_test():
    x = np.linspace(0,5,100)
    nb_bases = 25
    basis = create_spline_basis(x, nb_bases)
    assert_equals(basis.condition_number < 50, True)


def x_ref_test():
    x = np.linspace(0,5,100)
    x_ref = np.linspace(1,4,100)
    nb_bases = 25
    basis = create_spline_basis(x, nb_bases, x_ref=x_ref, \
            reg_coefs=[1e-3, 1e-3, 1e-3])
    assert_equals(len(basis.valid_idx), 60)
    assert_equals(basis.B_.shape, (205, 25))


def augment_test():
    x = np.linspace(0,5,100)
    nb_bases = 25
    basis = create_spline_basis(x, nb_bases, reg_coefs=[1e-3, 1e-3, 1e-3])

    y = np.random.rand(100)
    y_aug = basis.augment(y)

    assert_equals(y_aug.shape, (325,))
    assert_array_almost_equal_nulp(y_aug, np.r_[y, np.zeros(225)])


def brick_size_test_01():
    x = np.linspace(0,5,100)
    nb_bases = 25
    basis = create_spline_basis(x, nb_bases, reg_coefs=[1e-3, 0, 0])

    y = np.random.rand(100)
    y_aug = basis.augment(y)

    assert_equals(basis.B_.shape, (100+25, 25))
    assert_equals(y_aug.shape, (100+25,))


def brick_size_test_02():
    x = np.linspace(0,5,100)
    nb_bases = 25
    basis = create_spline_basis(x, nb_bases, reg_coefs=[0, 1, 0])

    y = np.random.rand(100)
    y_aug = basis.augment(y)

    assert_equals(basis.B_.shape, (100+100, 25))
    assert_equals(y_aug.shape, (100+100,))


def brick_size_test_03():
    x = np.linspace(0,5,100)
    nb_bases = 25
    basis = create_spline_basis(x, nb_bases, reg_coefs=[0, 0, 1])

    y = np.random.rand(100)
    y_aug = basis.augment(y)

    assert_equals(basis.B_.shape, (100+100, 25))
    assert_equals(y_aug.shape, (100+100,))


def brick_size_test_04():
    x = np.linspace(0,5,100)
    nb_bases = 25
    basis = create_spline_basis(x, nb_bases, reg_coefs=[1, 1, 1])

    y = np.random.rand(100)
    y_aug = basis.augment(y)

    assert_equals(basis.B_.shape, (100+225, 25))
    assert_equals(y_aug.shape, (100+225,))
