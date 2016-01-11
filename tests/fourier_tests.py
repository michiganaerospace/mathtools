'''Test basis generation algorithms for Fourier bases.'''
from mathtools.fit import *
import mathtools.refs.FS_Tools as fs
from mathtools.utils import *
from mathtools.fourier import *
from nose.tools import assert_equals, assert_almost_equals, assert_raises
from numpy.testing import assert_array_almost_equal_nulp, assert_array_equal
from nose import with_setup
import numpy as np


# BEGIN TESTS ------------------------------------------------------

# Can you create a Fourier basis? Let's find out.
def create_a_fourier_basis_test():
    x = np.linspace(0,5,100)
    B = fourier_basis(x, 25)
    assert_equals(B.shape, (100, 25))


def basis_matches_fs_tools_implementation_test():
    x = np.linspace(0,5,100)
    x_ = map_to_interval(x, [0,1])
    nb_bases = 51
    nb_terms = 25
    B = fourier_basis(x, nb_bases)
    B_fs = fs.FS_basis(x_, 1.0, nb_terms)
    assert_array_almost_equal_nulp(B, B_fs)


def create_a_derivative_fourier_basis_test():
    x = np.linspace(0,5,100)
    B = d_fourier_basis(x, 25)
    assert_equals(B.shape, (100, 25))


def d_basis_matches_fs_tools_implementation_test():
    x = np.linspace(0,5,100)
    x_ = map_to_interval(x, [0,1])
    nb_bases = 51
    nb_terms = 25
    dB = d_fourier_basis(x, nb_bases)
    dB_fs = fs.FS_dbasis(x_, 1.0, nb_terms)
    assert_array_almost_equal_nulp(dB, dB_fs)


def create_a_second_derivative_fourier_basis_test():
    x = np.linspace(0,5,100)
    d2B = d_fourier_basis(x, 25)
    assert_equals(d2B.shape, (100, 25))


def d2_basis_matches_fs_tools_implementation_test():
    x = np.linspace(0,5,100)
    x_ = map_to_interval(x, [0,1])
    nb_bases = 51
    nb_terms = 25
    d2B = d2_fourier_basis(x, nb_bases)
    d2B_fs = fs.FS_d2basis(x_, 1.0, nb_terms)
    assert_array_almost_equal_nulp(d2B, d2B_fs)


def create_fourier_test_01():
    x = np.linspace(0,5,100)
    nb_bases = 25
    basis = create_fourier_basis(x, nb_bases)
    assert_equals(basis.nb_bases, 25)


def create_fourier_shape_test():
    x = np.linspace(0,5,100)
    nb_bases = 25
    basis = create_fourier_basis(x, nb_bases)
    assert_equals(basis.B.shape, (100, 25))
    assert_equals(basis.dB.shape, (100, 25))
    assert_equals(basis.d2B.shape, (100, 25))
    assert_equals(basis.B_.shape, (325, 25))


def create_fourier_shape_test():
    x = np.linspace(0,5,100)
    nb_bases = 25
    basis = create_fourier_basis(x, nb_bases, reg_coefs=[1e-3, 1e-3, 1e-3])
    assert_equals(basis.inverse.shape, (25, 325))


def create_fourier_condition_number_test():
    x = np.linspace(0,5,100)
    nb_bases = 25
    basis = create_fourier_basis(x, nb_bases)
    assert_equals(basis.condition_number < 50, True)


def x_ref_test():
    x = np.linspace(0,5,100)
    x_ref = np.linspace(1,4,100)
    nb_bases = 25
    basis = create_fourier_basis(x, nb_bases, x_ref=x_ref, \
            reg_coefs=[1e-3, 1e-3, 1e-3])
    assert_equals(len(basis.valid_idx), 60)
    assert_equals(basis.B_.shape, (205, 25))


def augment_test():
    x = np.linspace(0,5,100)
    nb_bases = 25
    basis = create_fourier_basis(x, nb_bases, reg_coefs=[1e-3, 1e-3, 1e-3])

    y = np.random.rand(100)
    y_aug = basis.augment(y)

    assert_equals(y_aug.shape, (325,))
    assert_array_almost_equal_nulp(y_aug, np.r_[y, np.zeros(225)])


def brick_size_test_01():
    x = np.linspace(0,5,100)
    nb_bases = 25
    basis = create_fourier_basis(x, nb_bases, reg_coefs=[1e-3, 0, 0])

    y = np.random.rand(100)
    y_aug = basis.augment(y)

    assert_equals(basis.B_.shape, (100+25, 25))
    assert_equals(y_aug.shape, (100+25,))


def brick_size_test_02():
    x = np.linspace(0,5,100)
    nb_bases = 25
    basis = create_fourier_basis(x, nb_bases, reg_coefs=[0, 1, 0])

    y = np.random.rand(100)
    y_aug = basis.augment(y)

    assert_equals(basis.B_.shape, (100+100, 25))
    assert_equals(y_aug.shape, (100+100,))


def brick_size_test_03():
    x = np.linspace(0,5,100)
    nb_bases = 25
    basis = create_fourier_basis(x, nb_bases, reg_coefs=[0, 0, 1])

    y = np.random.rand(100)
    y_aug = basis.augment(y)

    assert_equals(basis.B_.shape, (100+100, 25))
    assert_equals(y_aug.shape, (100+100,))


def brick_size_test_04():
    x = np.linspace(0,5,100)
    nb_bases = 25
    basis = create_fourier_basis(x, nb_bases, reg_coefs=[1, 1, 1])

    y = np.random.rand(100)
    y_aug = basis.augment(y)

    assert_equals(basis.B_.shape, (100+225, 25))
    assert_equals(y_aug.shape, (100+225,))


if __name__ == '__main__':
    
    x = np.linspace(0,5,100)
    x_ref = np.linspace(1,4,100)
    nb_bases = 25
    basis = create_fourier_basis(x, nb_bases, x_ref=x_ref, \
            reg_coefs=[1e-3, 1e-3, 1e-3])
    assert_equals(len(basis.valid_idx), 60)
    assert_equals(basis.B_.shape, (205, 25))
