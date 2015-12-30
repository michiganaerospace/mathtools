'''Create bases for useful functions, including Fourier series, Legendre 
polynomials, and cubic splines.'''
import numpy as np
from mathtools.utils import map_to_interval


# LEGENDRE Polynomials ----------------------------------------------------

def legendre_basis(x, nb_bases):
    '''Define the Legendre polynomial basis.
    INPUTS
        x - array_like
            The samples on which we will sample the basis vectors. The vectors
            x must live on the interval [-1, 1]. In all cases, x will be mapped
            into this interval.
        nb_bases - int
            The number of basis vectors to create.
    OUTPUTS
        B - array_like
            A len(x) by nb_bases array. The columns of the array correspond to 
            the Legendre polynomial basis vectors. 
    '''
    # A little error checking.
    if nb_bases < 2:
        raise ValueError('Number of bases must satisfy nb_bases >= 2')

    # Map initial x --> x_ into the interval [-1,1]
    x_ = map_to_interval(x, [-1,1])
    nb_samples = len(x_)

    # Initialize the basis matrix, etc.
    B = np.zeros((nb_samples, nb_bases), dtype='double')
    p_km1 = np.ones(nb_samples, dtype='double') 
    p_k = x_

    # First two bases are ready to go.
    B[:,0] = p_km1
    B[:,1] = p_k

    # Recursively compute the remaining basis functions.
    for k in np.r_[2:nb_bases]:
        m = k - 1
        c = (2.0 * m + 1.0)
        p_kp1 = (c * x_ * p_k - m * p_km1)/np.double(m+1)
        B[:,k] = p_kp1
        p_km1 = p_k
        p_k = p_kp1

    return B


def d_legendre_basis(x, nb_bases):
    '''Define the first derivative of the legendre basis.
    INPUTS
        x - array_like
            The samples on which we will sample the basis vectors. The vectors
            x must live on the interval [-1, 1]. In all cases, x will be mapped
            into this interval.
        nb_bases - int
            The number of basis vectors to create.
    OUTPUTS
        B - array_like
            A len(x) by nb_bases array. The columns of the array correspond to 
            the the first derivative of the Legendre polynomial basis vectors. 
    '''
    # A little error checking.
    if nb_bases < 2:
        raise ValueError('Number of bases must satisfy nb_bases >= 2')

    # Map initial x --> x_ into the interval [-1,1]
    x_ = map_to_interval(x, [-1,1])
    nb_samples = len(x_)

    # Initialize the basis matrix.
    dB      = np.zeros((nb_samples, nb_bases), dtype='double')
    p_km1   = np.ones(nb_samples, dtype='double')
    p_k     = x_
    dp_km1  = np.zeros(nb_samples, dtype='double')
    dp_k    = np.ones(nb_samples, dtype='double')

    # First two basis vectors for free!
    dB[:,0] = dp_km1
    dB[:,1] = dp_k
    
    # Define the remaining basis vectors using recursion.
    for k in np.r_[2:nb_bases]:
        m       = k - 1
        c       = (2.0 * m + 1.0)
        p_kp1   = (c * x_ * p_k - m * p_km1)/np.double(m+1)
        dp_kp1  = (c * (p_k + x_*dp_k)- m*dp_km1)/np.double(m+1)
        dB[:,k] = dp_kp1
        p_km1   = p_k
        p_k     = p_kp1
        dp_km1  = dp_k
        dp_k    = dp_kp1

    return dB



def d2_legendre_basis(x, nb_bases):
    '''Define the second derivative of the legendre basis.
    INPUTS
        x - array_like
            The samples on which we will sample the basis vectors. The vectors
            x must live on the interval [-1, 1]. In all cases, x will be mapped
            into this interval.
        nb_bases - int
            The number of basis vectors to create.
    OUTPUTS
        B - array_like
            A len(x) by nb_bases array. The columns of the array correspond to 
            the the second derivative of the Legendre polynomial basis vectors. 
    '''
    # A little error checking.
    if nb_bases < 2:
        raise ValueError('Number of bases must satisfy nb_bases >= 2')

    # Map initial x --> x_ into the interval [-1,1]
    x_ = map_to_interval(x, [-1,1])
    nb_samples = len(x_)

    # Initialize the basis matrix.
    d2B = np.zeros((nb_samples, nb_bases), dtype='double')
    p_km1 = np.ones(nb_samples, dtype='double')
    p_k = x_
    dp_km1 = np.zeros(nb_samples, dtype='double')
    dp_k = np.ones(nb_samples, dtype='double')



    # d2Lbasis = zeros((N,Nterms),dtype='double')
    # Pkm1 = ones(N,dtype='double')
    # Pk = x
    # dPkm1 = zeros(N,dtype='double')
    # dPk = ones(N,dtype='double')

    # d2Pkm1 = zeros(N,dtype='double')
    # d2Pk = zeros(N,dtype='double')

    
    # d2Lbasis[:,1] = d2Pkm1
    # d2Lbasis[:,2] = d2Pk
    # for k in r_[2:Nterms]:
    #     m = k-1
    #     c = (2.0*m+1.0)
    #     Pkp1 = (c*x*Pk - m*Pkm1)/double(m+1)
    #     dPkp1 = (c*(Pk+x*dPk) - m*dPkm1)/double(m+1)
    #     d2Pkp1 = (c*(dPk + x*d2Pk + dPk) - m*d2Pkm1)/double(m+1)

    #     d2Lbasis[:,k] = d2Pkp1

    #     Pkm1 = Pk
    #     Pk = Pkp1
    #     dPkm1 = dPk
    #     dPk = dPkp1        
    #     d2Pkm1 = d2Pk
    #     d2Pk = d2Pkp1        
