'''Implements cubic spline basis functions.'''
import numpy as np
from mathtools.utils import Struct, pseudoinverse


def cubic_spline_basis_unit_interval(x):
    '''Compute four basis vectors over unit interval for x for a cubic spline.
    ARGUMENTS
        x - array_like
            The domain on which we wish to fit the data.
    OUTPUTS
        B - array_like 
            The four basis vectors defining the cubic spline on the unit 
            interval -- an (nb_samples x 4) array.
    '''
    # Initialize basis.
    nb_samples = len(x)
    B = np.zeros((nb_samples, 4), dtype='double')

    # Set up variables.
    u = 1.0 - x
    v = x
    c = 1.0/6.0
    q = 1.0 + (u*v)

    # Fill up the basis matrix!
    B[:,0] = c * (u**3) 
    B[:,1] = c + 0.5 * u * q
    B[:,2] = c + 0.5 * v * q
    B[:,3] = c * (v**3)

    return B


def d_cubic_spline_basis_unit_interval(x):
    '''Compute four basis vectors over unit interval for x for the derivative
       of the cubic spline.
    ARGUMENTS
        x - array_like
            The domain on which we wish to fit the data.
    OUTPUTS
        dB - array_like 
            The four basis vectors defining the cubic spline on the unit 
            interval -- an (nb_samples x 4) array.
    '''
    # Initialize basis.
    nb_samples = len(x)
    dB = np.zeros((nb_samples, 4), dtype='double')

    # Set up variables.
    u = 1.0 - x
    v = x
    dudx = -1.0
    dvdx = 1.0
    q = 1.0 + (u*v)

    dqdu = v
    dqdv = u

    # Fill up the basis matrix!
    dB[:,0] = 0.5 * (u**2) * dudx 
    dB[:,1] = 0.5 * ((u*dqdu + q)* dudx + (u*dqdv)*dvdx)
    dB[:,2] = 0.5 * ((v*dqdu)*dudx + (v*dqdv + q)*dvdx)
    dB[:,3] = 0.5 * (v**2)*dvdx

    return dB


def d2_cubic_spline_basis_unit_interval(x):
    '''Compute four basis vectors over unit interval for x for the second
       derivative of the cubic spline.
    ARGUMENTS
        x - array_like
            The domain on which we wish to fit the data.
    OUTPUTS
        d2B - array_like 
            The four basis vectors defining the cubic spline on the unit 
            interval -- an (nb_samples x 4) array.
    '''
    # Initialize basis.
    nb_samples = len(x)
    d2B = np.zeros((nb_samples, 4), dtype='double')

    # Set up variables.
    u = 1.0 - x
    v = x
    dudx = -1.0
    dvdx = 1.0
    q = 1.0 + (u*v)

    dqdx = dudx*v + u * dvdx
    dqdu = v
    dqdv = u
    d2qdudx = dvdx
    d2qdvdx = dudx

    # Fill up the basis matrix!
    d2B[:,0] = u*(dudx**2) 
    d2B[:,1] = 0.5*((dudx*dqdu + u*d2qdudx + dqdx)*dudx + \
            (dudx*dqdv+u*d2qdvdx)*dvdx)
    d2B[:,2] = 0.5*((dvdx*dqdu+v*d2qdudx)*dudx + \
            (dvdx*dqdv + v*d2qdvdx + dqdx)*dvdx)
    d2B[:,3] = v*(dvdx**2)

    return d2B


def uniform_knots(x, nb_knots):
    '''Compute four basis vectors over unit interval for x for the second
       derivative of the cubic spline.
    ARGUMENTS
        x - array_like
            The domain across which we will distribute the knots.
    OUTPUTS
        knots - array_like    
            An nb_knots length array of knot positions.
    '''

    # Distribute knots uniformly across interval spanned by x. NOte that we 
    # are dividing by (nb_knots - 3), which ensures the knots fully contain
    # the support of x.
    knots = x.min() + (x.max() - x.min()) * \
            np.r_[-1:nb_knots-1]/np.double(nb_knots-3)
    return knots


def cubic_spline_basis(x, nb_knots):
    '''Compute the cubic spline basis for domain x using nb_knots uniformly
       distributed knots.
    ARGUMENTS
        x - array_like
            The domain on which we wish to create the basis; an nb_samples 
            length array.
        nb_knots - int
            The number of uniformly distributed knots to use in constructing
            the basis.
    OUTPUTS
        B - array_like 
            The nb_samples x nb_knots array containing (column) cubic spline
            basis vectors.
    '''
    # Define variables.
    nb_samples = len(x)
    knots = uniform_knots(x, nb_knots)

    # Initialize basis.
    B = np.zeros((nb_samples, nb_knots), dtype='double')

    for k in np.r_[1:nb_knots-2]:
        idx_k = np.nonzero( (x>= knots[k]) & (x<=knots[k+1]) )

        # Map region between knots to the unit interval.
        t = (x[idx_k] - knots[k])/(knots[k+1]-knots[k])

        # Find the local cubic spline basis.
        local_B = cubic_spline_basis_unit_interval(t)

        # Fill appropriate columns of full basis.
        B[idx_k, k-1]   = local_B[:,0]
        B[idx_k, k]     = local_B[:,1]
        B[idx_k, k+1]   = local_B[:,2]
        B[idx_k, k+2]   = local_B[:,3]

    return B


def d_cubic_spline_basis(x, nb_knots):
    '''Compute the cubic spline basis for domain x using nb_knots uniformly
       distributed knots.
    ARGUMENTS
        x - array_like
            The domain on which we wish to create the basis; an nb_samples 
            length array.
        nb_knots - int
            The number of uniformly distributed knots to use in constructing
            the basis.
    OUTPUTS
        dB - array_like 
            The nb_samples x nb_knots array containing (column) derivative
            cubic spline basis vectors.
    '''
    # Define variables.
    nb_samples = len(x)
    knots = uniform_knots(x, nb_knots)

    # Initialize basis.
    dB = np.zeros((nb_samples, nb_knots), dtype='double')

    for k in np.r_[1:nb_knots-2]:
        idx_k = np.nonzero( (x>= knots[k]) & (x<=knots[k+1]) )

        # Map region between knots to the unit interval.
        h = knots[k+1] - knots[k]
        t = (x[idx_k] - knots[k])/h

        # Find the local cubic spline basis.
        local_dB = d_cubic_spline_basis_unit_interval(t)

        # Fill appropriate columns of full basis.
        dB[idx_k, k-1]   = local_dB[:,0]/h
        dB[idx_k, k]     = local_dB[:,1]/h
        dB[idx_k, k+1]   = local_dB[:,2]/h
        dB[idx_k, k+2]   = local_dB[:,3]/h

    return dB


def d2_cubic_spline_basis(x, nb_knots):
    '''Compute the cubic spline basis for domain x using nb_knots uniformly
       distributed knots.
    ARGUMENTS
        x - array_like
            The domain on which we wish to create the basis; an nb_samples 
            length array.
        nb_knots - int
            The number of uniformly distributed knots to use in constructing
            the basis.
    OUTPUTS
        d2B - array_like 
            The nb_samples x nb_knots array containing (column) second
            derivative cubic spline basis vectors.
    '''
    # Define variables.
    nb_samples = len(x)
    knots = uniform_knots(x, nb_knots)

    # Initialize basis.
    d2B = np.zeros((nb_samples, nb_knots), dtype='double')

    for k in np.r_[1:nb_knots-2]:
        idx_k = np.nonzero( (x>= knots[k]) & (x<=knots[k+1]) )

        # Map region between knots to the unit interval.
        h = knots[k+1] - knots[k]
        h2 = h**2
        t = (x[idx_k] - knots[k])/h

        # Find the local cubic spline basis.
        local_d2B = d2_cubic_spline_basis_unit_interval(t)

        # Fill appropriate columns of full basis.
        d2B[idx_k, k-1]   = local_d2B[:,0]/h2
        d2B[idx_k, k]     = local_d2B[:,1]/h2
        d2B[idx_k, k+1]   = local_d2B[:,2]/h2
        d2B[idx_k, k+2]   = local_d2B[:,3]/h2

    return d2B


def create_spline_basis(x, nb_bases, reg_coefs=[0,0,0], x_ref=None):
    '''Build a spline basis object.
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
            A reference domain. This is useful for resampling data. It ensures
            we avoid attempting to fit data outside of the original domain.
    OUTPUTS
        basis - Struct object
            A struct object containing the following methods and properties:
                - nb_bases: The number of basis vectors.
                - reg_coefs: A list of regularization coefficients.
                - x: The domain over which the basis will be defined.
                - valid_idx: indices of the domain that are valid -- that is,
                  those values inside the reference domain, if specified. If
                  not specified, the entire vector is valid.
                - B: Spline basis vectors (as columns).
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

    # Is there a reference domain? Note in the case of splines, there is no
    # need to scale data to live in the unit interval.
    if (x_ref is not None):
        x_min = x_ref.min()
        x_max = x_ref.max()
        basis.valid_idx = np.nonzero((x_ >= x_min) * (x_ <= x_max))[0]
        x_ = x_[basis.valid_idx]

    # Build bases and the 'brick'.
    basis.B      = cubic_spline_basis(x_, nb_bases)
    I            = np.eye(nb_bases)
    basis.dB     = d_cubic_spline_basis(x_, nb_bases)
    basis.d2B    = d2_cubic_spline_basis(x_, nb_bases)

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
