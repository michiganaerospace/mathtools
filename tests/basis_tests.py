from mathtools.fit import *
import mathtools.FS_Tools as fs
from mathtools.utils import *
from mathtools.bases import *
from nose.tools import assert_equals, assert_almost_equals, assert_raises
from numpy.testing import assert_array_almost_equal_nulp, assert_array_equal
from nose import with_setup
import numpy as np


# BEGIN TESTS ------------------------------------------------------

# Can you create a Legendre basis? Let's find out.
def create_a_legendre_basis():
    x = np.linspace(0,5,100)
    B = legendre_basis(x, 25)
    assert_equals(B.shape, (100, 25))


def matches_fs_tools_implementation_test():
    x = np.linspace(0,5,100)
    nb_bases = 25
    B = legendre_basis(x, nb_bases)
    B_fs = fs.Legendre_basis(nb_bases, x)
    print (B-B_fs).max()
    assert((B - B_fs).max() < 1e-13)
    assert_array_almost_equal_nulp(B, B_fs, 40000)


def create_derivative_legendre_basis():
    x= linspace(0,5,100)
    dB = d_legendre_basis(x, 25)
    assert_equals(dB.shape, (100, 25))


def dl_matches_fs_tools_test():
    x = np.linspace(0,5,100)
    nb_bases = 25
    dB = d_legendre_basis(x, nb_bases)
    dB_fs = fs.dLegendre_basis(nb_bases, x)
    assert_array_almost_equal_nulp(dB, dB_fs, 40000)
    

if __name__ == '__main__':
    x = np.linspace(0,5,100)
    nb_bases = 25
    B = legendre_basis(x, nb_bases)
    B_fs = fs.Legendre_basis(nb_bases, x)
