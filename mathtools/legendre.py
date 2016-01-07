'''Create bases for useful functions, including Fourier series, Legendre 
polynomials, and cubic splines.'''
import numpy as np
from mathtools.utils import map_to_interval, pseudoinverse, Struct


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
        d2B - array_like
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
    d2B         = np.zeros((nb_samples, nb_bases), dtype='double')
    p_km1       = np.ones(nb_samples, dtype='double')
    p_k         = x_
    dp_km1      = np.zeros(nb_samples, dtype='double')
    dp_k        = np.ones(nb_samples, dtype='double')
    d2p_km1     = np.zeros(nb_samples, dtype='double')
    d2p_k       = np.zeros(nb_samples, dtype='double')

    # Add the first two bases.
    d2B[:,1]    = d2p_km1
    d2B[:,2]    = d2p_k

    # Define remaining basis vectors via recursion.
    for k in np.r_[2:nb_bases]:
        m       = k-1
        c       = (2.0 * m + 1.0)
        p_kp1   = (c * x_ * p_k - m * p_km1)/np.double(m+1)
        dp_kp1  = (c * (p_k + x_ * dp_k) - m * dp_km1)/np.double(m+1)
        d2p_kp1 = (c * (dp_k + x_ * d2p_k + dp_k) - m * d2p_km1)/np.double(m+1)

        d2B[:,k]    = d2p_kp1
        p_km1       = p_k
        p_k         = p_kp1
        dp_km1      = dp_k
        dp_k        = dp_kp1
        d2p_km1     = d2p_k
        d2p_k       = d2p_kp1
    
    return d2B


def create_legendre_basis(x, nb_bases, reg_coefs=[0,0,0]):
    '''Build legendre polynomial bases.
    INPUTS
        x - array_like
            An array of points -- the domain on which we will build the basis.
        nb_bases - int
            The number of basis vectors to generate.
        reg_coefs - array_like (default: [0.0, 0.0, 0.0])
            An array or list of three numerical coefficients that specify 
            the regularization penalty for the magnitude of coefficients, as
            well as the the magnitude of the first and second derivatives.
    OUTPUTS
        basis - Struct object
            A struct object containing the following fields:
                - nb_bases: the number of basis vectors
                - reg_coefs: a list of regularization coefficients
                - x: the domain over which the basis is defined
                - B: Legendre basis vectors (as columns)
                - dB: derivative of basis vectors in B (times reg_coefs[1])
                - d2B: second derivative of basis vectors in B (times 
                       reg_coefs[2])
                - I: identity matrix (times reg_coefs[0])
                - B_: the 'brick', a concatenation of B, I, dB, and d2B
                - inverse: the pseudo-inverse of the brick, B_
                - condition_number: the condition number of the inverse of the
                  brick.
    '''
    # Build a structure to hold the data.
    basis           = Struct()
    basis.nb_bases  = nb_bases 
    basis.reg_coefs = reg_coefs
    basis.x         = x
    
    # Build bases and the 'brick'.
    basis.B      = legendre_basis(x, nb_bases)
    basis.I      = reg_coefs[0] * np.eye(nb_bases)
    basis.dB     = reg_coefs[1] * d_legendre_basis(x, nb_bases)
    basis.d2B    = reg_coefs[2] * d2_legendre_basis(x, nb_bases)

    # Create the 'brick' by stacking these bases on top of one another.
    basis.B_     = np.r_[basis.B, basis.I, basis.dB, basis.d2B] 

    # Find the inverse of the brick. Keep the condition number around, too.
    basis.inverse, basis.condition_number = pseudoinverse(basis.B_, True)

    return basis
