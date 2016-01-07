import numpy as np
from mathtools.fit import *
from mathtools.utils import *


if __name__ == '__main__':
    
    # Generate some matrices.
    M = np.random.rand(15,10)
    N = np.random.rand(10,15)

    # Compute the pseudoinverse using the numpy's function.
    Minv = np.linalg.pinv(M)
    Ninv = np.linalg.pinv(N)

    # Use our version of the pseudoinverse.
    M_inv_mt, condition_M = pseudoinverse(M, True)
    N_inv_mt, condition_N = pseudoinverse(N, True)



