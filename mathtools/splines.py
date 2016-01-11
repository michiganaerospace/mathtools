import numpy as np


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

