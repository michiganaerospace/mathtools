'''Utility functions for mathtools... tools.'''
import numpy as np


def validate_type(var, white_list, var_name='Input'):
    '''Ensure that the variable type is the supplied whitelist. 
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
   

def scale_to_interval(x, interval):
    '''Shift and scale vector so that its elements live in the specified 
       interval.
    INPUTS
        x - array_like
            The input vector.
        interval - array_like
            The interval, [a,b], into which we will shift and stretch the
            elements of vector <x>.
    OUTPUTS
        x_ - array_like
            A new vector, with the same dimensions as <x>, but each of whose 
            elements now live the within the interval specified in <interval>.
    '''

    # Ensure we are working with numpy arrays or lists.
    validate_type(x, [list, np.ndarray])
    validate_type(interval, [list, np.ndarray])

    if type(interval) == list:
        interval = np.array(interval)
    
    # Start with the original vector.
    x_ = np.array(x)

    # Apply appropriate shift.
    shift = x.mean() - interval.mean()
    x_ -= shift

    return x_
    

