import numpy as np
from mathtools.utils import map_to_interval


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
        dB[:,k+1] = -2.0 * np.sin(2.0 * np.pi * (k+1) * freq * x)
        dB[:,k+1+nb_terms] = 2.0 * np.cos(2.0 * np.pi * (k+1) * freq * x)

    return dB


def d2_fourier_basis(x, nb_bases, freq=1.0):
    '''Construct the second derivative of the Fourier series basis with phase
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


def FS_dbasis(R2,alpha,NFterms):
    N = len(R2)
    dbasis = zeros((N,2*NFterms+1),dtype='double')
    dbasis[:,0] = 0.0
    for k in r_[0:NFterms]:
        dbasis[:,k+1] = -2.0*sin(2.0*pi*(k+1)*alpha*R2)*(k+1)
        dbasis[:,k+1+NFterms] = 2.0*cos(2.0*pi*(k+1)*alpha*R2)*(k+1)

    return dbasis
