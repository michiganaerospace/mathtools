from mathtools.fit import *
from mathtools.utils import *
from mathtools.legendre import *
from nose.tools import assert_equals, assert_almost_equals, assert_raises
from numpy.testing import assert_array_almost_equal_nulp, assert_array_equal
from nose import with_setup
import os, glob
import numpy as np


# BEGIN TESTS ------------------------------------------------------

def validate_type_test():
    valid_types = [list]
    is_good = validate_type([1,2,3], valid_types)
    assert_equals(is_good, True)


def fail_single_type_validation_test():
    valid_types = [list]
    assert_raises(ValueError, validate_type, 42, valid_types)


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


def vessel_init_test():
    v = Vessel()
    v.name = 'Matthew Lewis'
    assert_equals(v.name, 'Matthew Lewis')


def vessel_filename_test():
    v = Vessel('myfilename.dat')
    assert_equals(v.current_filename, 'myfilename.dat')


def vessel_save_test():
    v = Vessel('myfilename.dat')
    v.payload = [1,2,3]
    v.save()
    assert_equals(glob.glob('*.dat')[0], 'myfilename.dat')

    # Clean up after ourselves.
    for filename in glob.glob('*.dat'):
        os.remove(filename)


def vessel_save_error_test():
    v = Vessel() 
    v.payload = [3,1,4,1,5]
    assert_raises(ValueError, v.save)


def vessel_load_test():
    v = Vessel('myfilename.dat')
    v.payload = [1,2,3]
    v.save()

    g = Vessel(v.current_filename)
    assert_equals(g.payload, v.payload)

    # Clean up after ourselves.
    for filename in glob.glob('*.dat'):
        os.remove(filename)


def vessel_ingest_test():
    v = Vessel('myfilename.dat')
    data = {'name':'Matt', 'age': 39}
    v.ingest(data)
    assert_equals(v.name, 'Matt')
    assert_equals(v.age, 39)


def vessel_keys_test():
    v = Vessel()
    v.name = 'Matt'
    v.kids = ['Sophie', 'Claire', 'Hugo']
    v.ages = [11, 8, 0]
    assert_array_equal(v.keys, ['ages', 'kids', 'name', 'current_filename'])


def mahal_test():
    x = np.random.rand(100, 5)
    (d, mu, S) = mahal(x, return_stats=True)
    assert_equals(d.shape, (100,))
    assert_equals(mu.shape, (5,))
    assert_equals(S.shape, (5,5))


def mahal_simple_test():
    x = np.random.rand(100, 5)
    d = mahal(x)
    assert_equals(d.shape, (100,))


def mahal_with_mean_test():
    x = np.random.rand(100,10)
    mu = x.mean(0)
    S = np.cov(x.T)
    (d, mu_out, S_out) = mahal(x, mu, S, True)
    assert_array_almost_equal_nulp(mu, mu_out)
    assert_array_almost_equal_nulp(S, S_out)

    

