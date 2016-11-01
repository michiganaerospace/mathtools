'''Test basis generation algorithms for cubic splines, Legendre polynomials,
and Fourier series.'''
from mathtools.fit import *
import mathtools.refs.FS_Tools as fs
from mathtools.utils import *
from mathtools.legendre import *
from nose.tools import assert_equals, assert_almost_equals, assert_raises
from numpy.testing import assert_array_almost_equal_nulp, assert_array_equal
from nose import with_setup
import numpy as np


# BEGIN TESTS ------------------------------------------------------

# Can you create a Legendre basis? Let's find out.
def create_a_legendre_basis_test():
    x = np.linspace(0,5,100)
    B = legendre_basis(x, 25)
    assert_equals(B.shape, (100, 25))


def matches_fs_tools_implementation_test():
    x = np.linspace(0,5,100)
    nb_bases = 25
    B = legendre_basis(x, nb_bases)
    B_fs = fs.Legendre_basis(nb_bases, x)
    # print (B-B_fs).max()
    assert((B - B_fs).max() < 1e-13)
    assert_array_almost_equal_nulp(B, B_fs)


def create_derivative_legendre_basis_test():
    x = np.linspace(0,5,100)
    dB = d_legendre_basis(x, 25)
    assert_equals(dB.shape, (100, 25))


def dl_matches_fs_tools_test():
    x = np.linspace(0,5,100)
    nb_bases = 25
    dB = d_legendre_basis(x, nb_bases)
    dB_fs = fs.dLegendre_basis(nb_bases, x)
    assert_array_almost_equal_nulp(dB, dB_fs)
    

def create_second_derivative_legendre_basis():
    x = np.linspace(0,5,100)
    dB = d2_legendre_basis(x, 25)
    assert_equals(dB.shape, (100, 25))


def d2l_matches_fs_tools_test():
    x = np.linspace(0,5,100)
    nb_bases = 25
    d2B = d2_legendre_basis(x, nb_bases)
    d2B_fs = fs.d2Legendre_basis(nb_bases, x)
    assert_array_almost_equal_nulp(d2B, d2B_fs)


def throw_error_on_too_few_bases_test():
    x = np.linspace(0,5,100)
    nb_bases = 1
    assert_raises(ValueError, legendre_basis, x, nb_bases)
    assert_raises(ValueError, d_legendre_basis, x, nb_bases)
    assert_raises(ValueError, d2_legendre_basis, x, nb_bases)


def create_legendre_test_01():
    x = np.linspace(0,5,100)
    nb_bases = 25
    basis = create_legendre_basis(x, nb_bases)
    assert_equals(basis.nb_bases, 25)


def create_legendre_shape_test():
    x = np.linspace(0,5,100)
    nb_bases = 25
    basis = create_legendre_basis(x, nb_bases)
    assert_equals(basis.B.shape, (100, 25))
    assert_equals(basis.dB.shape, (100, 25))
    assert_equals(basis.d2B.shape, (100, 25))
    assert_equals(basis.B_.shape, (325, 25))


def create_legendre_shape_test():
    x = np.linspace(0,5,100)
    nb_bases = 25
    basis = create_legendre_basis(x, nb_bases, reg_coefs=[1e-3, 1e-3, 1e-3])
    assert_equals(basis.inverse.shape, (25, 325))


def create_legendre_condition_number_test():
    x = np.linspace(0,5,100)
    nb_bases = 25
    basis = create_legendre_basis(x, nb_bases)
    assert_equals(basis.condition_number < 50, True)


def x_ref_test():
    x = np.linspace(0,5,100)
    x_ref = np.linspace(1,4,100)
    nb_bases = 25
    basis = create_legendre_basis(x, nb_bases, x_ref=x_ref, \
            reg_coefs=[1e-3, 1e-3, 1e-3])
    assert_equals(len(basis.valid_idx), 60)
    assert_equals(basis.B_.shape, (205, 25))


def augment_test():
    x = np.linspace(0,5,100)
    nb_bases = 25
    basis = create_legendre_basis(x, nb_bases, reg_coefs=[1e-3, 1e-3, 1e-3])

    y = np.random.rand(100)
    y_aug = basis.augment(y)

    assert_equals(y_aug.shape, (325,))
    assert_array_almost_equal_nulp(y_aug, np.r_[y, np.zeros(225)])


def brick_size_test_01():
    x = np.linspace(0,5,100)
    nb_bases = 25
    basis = create_legendre_basis(x, nb_bases, reg_coefs=[1e-3, 0, 0])

    y = np.random.rand(100)
    y_aug = basis.augment(y)

    assert_equals(basis.B_.shape, (100+25, 25))
    assert_equals(y_aug.shape, (100+25,))


def brick_size_test_02():
    x = np.linspace(0,5,100)
    nb_bases = 25
    basis = create_legendre_basis(x, nb_bases, reg_coefs=[0, 1, 0])

    y = np.random.rand(100)
    y_aug = basis.augment(y)

    assert_equals(basis.B_.shape, (100+100, 25))
    assert_equals(y_aug.shape, (100+100,))


def brick_size_test_03():
    x = np.linspace(0,5,100)
    nb_bases = 25
    basis = create_legendre_basis(x, nb_bases, reg_coefs=[0, 0, 1])

    y = np.random.rand(100)
    y_aug = basis.augment(y)

    assert_equals(basis.B_.shape, (100+100, 25))
    assert_equals(y_aug.shape, (100+100,))


def brick_size_test_04():
    x = np.linspace(0,5,100)
    nb_bases = 25
    basis = create_legendre_basis(x, nb_bases, reg_coefs=[1, 1, 1])

    y = np.random.rand(100)
    y_aug = basis.augment(y)

    assert_equals(basis.B_.shape, (100+225, 25))
    assert_equals(y_aug.shape, (100+225,))
