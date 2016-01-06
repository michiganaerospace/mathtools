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

    # Compute the thing that does the stuff
    U, s, Vt = np.linalg.svd(N)
    V = Vt.T
    maxdim = len(s)
    M_inv = V[:,0:maxdim].dot(np.diag(1/s)).dot(U[:,0:maxdim].T)





