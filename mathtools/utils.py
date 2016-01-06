'''Utility functions for mathtools... tools.'''
import numpy as np


class Struct(object):
    '''Generic structure for conveniently holding properties.'''
    pass


def validate_type(var, white_list, var_name='Input'):
    '''Ensure that the variable type is in the supplied whitelist. 
    INPUTS
        var - variable or object
            A valid Python variable or object.
        white_list - list
            A list of allowable types. If type(<var>) is not in <white_list>,
            an error is thrown.
    OUTPUTS
        valid_type - boolean
            True if the variable is in the whitelist. Otherwise an error is 
            thrown with an appropriate error message.
    '''
    if type(var) in white_list:
        return True
    else:
        error_message = '%s must be of type: ' % var_name
        for idx, allowable_type in enumerate(white_list):
            if idx == len(white_list)-1:
                error_message += str(allowable_type)
            else:
                error_message += str(allowable_type) + ', or '
        raise ValueError(error_message)
   

def map_to_interval(x, interval, return_all=False):
    '''Shift and scale vector so that its elements live in the specified 
       interval.
    INPUTS
        x - array_like
            The input vector.
        interval - array_like
            The interval, [a,b], into which we will shift and stretch the
            elements of vector <x>.
        return_all - Boolean
            If True, the function returns scaled data, as well as the shift
            and scaling functions: 
                (x_, shift, scale)
            such that x_ = scale * (x - shift)
    OUTPUTS
        x_ - array_like
            A new vector, with the same dimensions as <x>, but each of whose 
            elements now live the within the interval specified in <interval>.
        shift - float
            The amount by which the data is shifted before it is scaled.
            This is returned only if return_all=True.
        scale - float
            The amount by which the data was scaled, after it was shifted. 
            This is returned only if return_all=True.
    '''

    # Ensure we are working with numpy arrays or lists.
    validate_type(x, [list, np.ndarray])
    validate_type(interval, [list, np.ndarray])

    if type(interval) == list:
        interval = np.array(interval)

    if len(interval) != 2:
        raise ValueError('The interval must be a size 2 vector: [a, b]')
    
    # Define the shift and scale parameters.
    scale = (interval.max() - interval.min())/(x.max() - x.min())
    shift = x.min() - interval.min()/scale

    # Scale & shift the vector.
    x_ = scale * (x - shift)

    if return_all:
        return x_, shift, scale
    else:
        return x_


def pseudoinverse(M):
    '''Find the pseudoinverse of the matrix M using singular value
       decomposition.
    INPUT
        M - array_like
            An (m x n) matrix whose pseudoinverse we want to find.
    OUTPUT
        pinv - array_like
            The Moore-Penrose pseudoinverse of the matrix M.
        cond - float
            The condition number of the inversion (i.e., the ratio of the
            largest to smallest singular values).
    '''
    # Compute the singular value decomposition.
    U, s, Vt = np.linalg.svd(N)
    V = Vt.T
    max_dim = len(s)
    M_pinv = V[:,:max_dim].dot(np.diag(1/s)).dot(U[:,:max_dim].T)
    condition_number = s.max()/s.min() 
