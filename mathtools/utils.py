'''Utility functions for mathtools tools...'''
import numpy as np


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

    # Ensure we are working with numpy arrays.
    if type(x) == list:
        x = np.array(x)

    if type(interval) == list:
        interval = np.array(interval)
    
    # Start with the original vector.
    x_ = np.array(x)

    # Apply appropriate shift.
    shift = x.mean() - interval.mean()
    x_ -= shift


    return x_
    

