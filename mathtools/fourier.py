import numpy as np
from mathtools.utils import map_to_interval, Struct, pseudoinverse


def fourier_basis(x, nb_bases, freq=1.0):
    '''Construct a Fourier series basis with phase 2*pi*freq*x.
    ARGUMENTS
        x - array_like
            The domain on which we are building the Fourier basis; an 
            nb_samples length array.
        nb_bases - int
            The number of Fourier bases to construct. Note that this is not the
            number of terms in the individual sin and cos series. This is the
            total number of basis vectors. nb_bases = 2 * nb_terms + 1.
        freq - float [default: 1.0]
            The fundamental frequency of the series.
    OUTPUT
        B - array_like
            An nb_samples x nb_bases size array whose columns are the Fourier
            basis vectors. 
    '''

    # Error checking -- nb_bases must be an odd number.
    if np.mod(nb_bases, 2) == 0:
        raise ValueError('Number of bases must be an odd number.')

    # Map the domain to the interval [0, 1].
    x_ = map_to_interval(x, [0, 1])

    # Initialize basis array.
    nb_terms = (nb_bases - 1)/2.0
    nb_samples = len(x_)
    B = np.zeros((nb_samples, 2*nb_terms+1), dtype='double')
    B[:,0] = 1.0

    # Construct the basis (column) vectors.
    for k in np.r_[0:nb_terms]:
        B[:,k+1] = 2.0*np.cos(2.0*np.pi*(k+1)*freq*x_)
        B[:,k+1+nb_terms] = 2.0*np.sin(2.0*np.pi*(k+1)*freq*x_)

    return B


def d_fourier_basis(x, nb_bases, freq=1.0):
    '''Construct the first derivative of the Fourier series basis with phase
       2*pi*freq*x.
    ARGUMENTS
        x - array_like
            The domain on which we are building the Fourier basis; an 
            nb_samples length array.
        nb_bases - int
            The number of Fourier bases to construct. Note that this is not the
            number of terms in the individual sin and cos series. This is the
            total number of basis vectors. nb_bases = 2 * nb_terms + 1.
        freq - float [default: 1.0]
            The fundamental frequency of the series.
    OUTPUT
        dB - array_like
            An nb_samples x nb_bases size array whose columns are the
            derivative of the Fourier basis vectors. 
    '''

    # Error checking -- nb_bases must be an odd number.
    if np.mod(nb_bases, 2) == 0:
        raise ValueError('Number of bases must be an odd number.')

    # Map the domain to the interval [0, 1].
    x_ = map_to_interval(x, [0, 1])

    # Initialize basis array.
    nb_terms = (nb_bases - 1)/2.0
    nb_samples = len(x_)
    dB = np.zeros((nb_samples, 2*nb_terms+1), dtype='double')
    dB[:,0] = 0.0

    # Construct the basis (column) vectors.
    for k in np.r_[0:nb_terms]:
        dB[:,k+1] = -2.0*np.sin(2.0*np.pi*(k+1)*freq*x_)*(k+1)
        dB[:,k+1+nb_terms] = 2.0*np.cos(2.0*np.pi*(k+1)*freq*x_)*(k+1)

    return dB


def d2_fourier_basis(x, nb_bases, freq=1.0):
    '''Construct the second derivative of the Fourier series basis with phase
       2*pi*freq*x.
    ARGUMENTS
        x - array_like
            The domain on which we are building the Fourier basis; an 
            nb_samples length array.
        nb_bases - int
            The number of Fourier bases to construct. NOTE that this is not the
            number of terms in the individual sin and cos series. This is the
            total number of basis vectors. nb_bases = 2 * nb_terms + 1.
        freq - float [default: 1.0]
            The fundamental frequency of the series.
    OUTPUT
        d2B - array_like
            An nb_samples x nb_bases size array whose columns are the
            derivative of the Fourier basis vectors. 
    '''

    # Error checking -- nb_bases must be an odd number.
    if np.mod(nb_bases, 2) == 0:
        raise ValueError('Number of bases must be an odd number.')

    # Map the domain to the interval [0, 1].
    x_ = map_to_interval(x, [0, 1])

    # Initialize basis array.
    nb_terms = (nb_bases - 1)/2.0
    nb_samples = len(x_)
    d2B = np.zeros((nb_samples, 2*nb_terms+1), dtype='double')
    d2B[:,0] = 0.0

    # Construct the basis (column) vectors.
    for k in np.r_[0:nb_terms]:
        d2B[:,k+1] = -2.0*np.cos(2.0*np.pi*(k+1)*freq*x_)*(k+1)**2
        d2B[:,k+1+nb_terms] = -2.0*np.sin(2.0*np.pi*(k+1)*freq*x_)*(k+1)**2

    return d2B


def create_fourier_basis(x, nb_bases, freq=1.0, reg_coefs=[0,0,0], x_ref=None):
    '''Build Fourier series basis object.
    INPUTS
        x - array_like
            An array of points -- the domain on which we will build the basis.
        nb_bases - int
            The number of basis vectors to generate.
        freq - float
            The fundamental frequency of the series.
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
            A struct object containing the following methods and properties:
                - nb_bases: The number of basis vectors.
                - reg_coefs: A list of regularization coefficients.
                - x: The domain over which the basis will be defined.
                - valid_idx: indices of the domain that are valid -- that is,
                  those values inside the reference domain, if specified. If
                  not specified, the entire vector is valid.
                - B: Legendre basis vectors (as columns).
                - dB: Derivative of basis vectors in B 
                - d2B: Second derivative of basis vectors in B 
                - B_: The 'brick', a concatenation of B, I, dB, and d2B.
                - inverse: The pseudo-inverse of the brick, B_.
                - condition_number: The condition number associated with the
                  brick inverse.
                - augment(y): A method that takes in an nb_samples length 
                  data vector, y, and returns a properly augmented data vector.
                  The vector is concatenated with the proper number of zeros
                  so that regularized least squares just works.
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
        x_ref_map, shift, scale = map_to_interval(x_ref, [0,1], \
                return_all=True)

        # Map current domain using reference domain scale and shift.
        x_ = scale * (x - shift)
        basis.valid_idx = np.nonzero((x_ >= 0) * (x_ <= 1))[0]
        x_ = x_[basis.valid_idx]

    # Build bases and the 'brick'.
    basis.B      = fourier_basis(x_, nb_bases, freq)
    I            = np.eye(nb_bases)
    basis.dB     = d_fourier_basis(x_, nb_bases, freq)
    basis.d2B    = d2_fourier_basis(x_, nb_bases, freq)

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
