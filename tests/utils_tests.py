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
