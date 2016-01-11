from mathtools.fit import *
from mathtools.utils import *
from mathtools.legendre import *
from nose.tools import assert_equals, assert_almost_equals, assert_raises
from numpy.testing import assert_array_almost_equal_nulp, assert_array_equal
from nose import with_setup
import numpy as np


# BEGIN TESTS ------------------------------------------------------

# Can we test for valid types?
def validate_type_test():
    valid_types = [list]
    is_good = validate_type([1,2,3], valid_types)
    assert_equals(is_good, True)


def fail_single_type_validation_test():
    valid_types = [list]
    assert_raises(ValueError, validate_type, 42, valid_types)


# Mapping data to a specified interval.
def scale_vector_data_to_unit_interval_test():
    x = np.random.rand(50)
    x_ = map_to_interval(x, [0, 1])
    assert_equals(x_.min(), 0)
    assert_almost_equals(x_.max(), 1)


def scale_to_straddling_test():
    x = np.random.rand(50)
    x_ = map_to_interval(x, [-3, 3])
    assert_almost_equals(x_.min(), -3)
    assert_almost_equals(x_.max(), 3)


def scale_to_negative_interval_test():
    x = np.random.rand(50)
    x_ = map_to_interval(x, [-3, -1])
    assert_almost_equals(x_.min(), -3)
    assert_almost_equals(x_.max(), -1)


def scale_to_positive_interval_test():
    x = np.random.rand(50)
    x_ = map_to_interval(x, [13, 21])
    assert_almost_equals(x_.min(), 13)
    assert_almost_equals(x_.max(), 21)


def return_all_test():
    x = np.random.rand(50)
    x_, shf, scl = map_to_interval(x, [0,1], return_all=True)
    assert_almost_equals(x_.min(), 0)
    assert_almost_equals(x_.max(), 1)


def pseudoinverse_test():
    np.random.seed(200)
    M = np.random.rand(15,10)
    M_inv = pseudoinverse(M)
    assert_array_almost_equal_nulp(M_inv, np.linalg.pinv(M), 500)
    assert_equals(M_inv.shape, (10,15))


def pseudoinverse_test_2():
    np.random.seed(200)
    M = np.random.rand(10,15)
    M_inv = pseudoinverse(M)
    assert_array_almost_equal_nulp(M_inv, np.linalg.pinv(M), 1500)
    assert_equals(M_inv.shape, (15,10))


def pseudoinverse_test_3():
    M = np.random.rand(10,15)
    M_inv, condition_number = pseudoinverse(M, True)
    assert_equals(condition_number<100, True)


