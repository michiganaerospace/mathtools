
'''Test basis generation algorithms for cubic splines, Legendre polynomials,
and Fourier series.'''
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
