'''Create bases for useful functions, including Fourier series, Legendre 
polynomials, and cubic splines.'''
import numpy as np
from mathtools.utils import map_to_interval, pseudoinverse, Struct
import pdb


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


def create_legendre_basis(x, nb_bases, reg_coefs=[0,0,0], x_ref=None):
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
        x_ref - array_like
            A reference domain. This is useful for resampling data.  It ensures
            that data is mapped to the interval [-1, 1] in the same way, and
            allows us to avoid attempting to fit data outside of the original
            domain.
    OUTPUTS
        basis - Struct object
            A struct object containing the following fields:
                - nb_bases: The number of basis vectors.
                - reg_coefs: A list of regularization coefficients.
                - x: The domain over which the basis will be defined.
                - valid_idx: indices of the domain that are valid -- that is,
                  those values inside the reference domain, if specified. If
                  not specified, the entire vector is valid.
                - B: Legendre basis vectors (as columns).
                - I: Identity matrix (times reg_coefs[0]).
                - dB: Derivative of basis vectors in B 
                - d2B: Second derivative of basis vectors in B 
                - B_: The 'brick', a concatenation of B, I, dB, and d2B.
                - inverse: The pseudo-inverse of the brick, B_.
                - condition_number: The condition number associated with the
                  brick inverse.
    '''
    # Build a structure to hold the data.
    basis           = Struct()
    basis.nb_bases  = nb_bases 
    basis.reg_coefs = reg_coefs
    basis.x         = x
    
    # Define the invalid indices (empty set to start).
    basis.valid_idx = np.arange(len(x)) 

    # This is default shift and scaling.
    x_ = x

    # Is there a reference domain? 
    if (x_ref is not None):
        x_ref_map, shift, scale = map_to_interval(x_ref, [-1,1], \
                return_all=True)

        # Map current domain using reference domain scale and shift.
        x_ = scale * (x - shift)
        basis.valid_idx = np.nonzero((x_ >= -1) * (x_ <= 1))[0]
        x_ = x_[basis.valid_idx]

    # Build bases and the 'brick'.
    basis.B      = legendre_basis(x_, nb_bases)
    I            = np.eye(nb_bases)
    basis.dB     = d_legendre_basis(x_, nb_bases)
    basis.d2B    = d2_legendre_basis(x_, nb_bases)

    # Create the 'brick' by stacking these bases on top of one another. Only
    # include those components that have nonzero regularization coefficients.
    # This can speed up the SVD computation.
    basis.B_ = basis.B
    if reg_coefs[0] > 0:
        basis.B_ = np.r_[basis.B_, reg_coefs[0] * I]
    if reg_coefs[1] > 0:
        basis.B_ = np.r_[basis.B_, reg_coefs[1] * basis.dB]
    if reg_coefs[2] > 0:
        basis.B_ = np.r_[basis.B_, reg_coefs[2] * basis.d2B]

    # Find the inverse of the brick. Keep the condition number around, too.
    basis.inverse, basis.condition_number = pseudoinverse(basis.B_, True)

    # Define the data augmention function.
    def augment(y):
        nb_zeros = basis.B_.shape[0] - len(y)
        return np.r_[y, np.zeros(nb_zeros)]

    # Attach augment function to the basis object.
    basis.augment = augment

    return basis
