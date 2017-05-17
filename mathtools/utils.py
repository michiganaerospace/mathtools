'''Utility functions for mathtools... tools.'''
import numpy as np
try: 
    import cPickle as pickle
except:
    import pickle
import glob


class Struct(object):
    '''Generic lightweight structure for conveniently holding properties.'''
    pass


def validate_type(var, white_list, var_name='Input'):
    '''Ensure that the variable type is in the supplied whitelist. 
    INPUTS
        var - variable or object
            A valid Python variable or object.
        white_list - list
            A list of allowable types. If type(<var>) is not in <white_list>,
            an error is thrown.
    OUTPUTS
        valid_type - boolean
            True if the variable is in the whitelist. Otherwise an error is 
            thrown with an appropriate error message.
    '''
    if type(var) in white_list:
        return True
    else:
        error_message = '%s must be of type: ' % var_name
        for idx, allowable_type in enumerate(white_list):
            if idx == len(white_list)-1:
                error_message += str(allowable_type)
            else:
                error_message += str(allowable_type) + ', or '
        raise ValueError(error_message)
   

def map_to_interval(x, interval, return_all=False):
    '''Shift and scale vector so that its elements live in the specified 
       interval.
    INPUTS
        x - array_like
            The input vector.
        interval - array_like
            The interval, [a,b], into which we will shift and stretch the
            elements of vector <x>.
        return_all - Boolean
            If True, the function returns scaled data, as well as the shift
            and scaling functions: 
                (x_, shift, scale)
            such that x_ = scale * (x - shift)
    OUTPUTS
        x_ - array_like
            A new vector, with the same dimensions as <x>, but each of whose 
            elements now live the within the interval specified in <interval>.
        shift - float
            The amount by which the data is shifted before it is scaled.
            This is returned only if return_all=True.
        scale - float
            The amount by which the data was scaled, after it was shifted. 
            This is returned only if return_all=True.
    '''

    # Ensure we are working with numpy arrays or lists.
    validate_type(x, [list, np.ndarray])
    validate_type(interval, [list, np.ndarray])

    if type(interval) == list:
        interval = np.array(interval)

    if len(interval) != 2:
        raise ValueError('The interval must be a size 2 vector: [a, b]')
    
    # Define the shift and scale parameters.
    scale = (interval.max() - interval.min())/(x.max() - x.min())
    shift = x.min() - interval.min()/scale

    # Scale & shift the vector.
    x_ = scale * (x - shift)

    if return_all:
        return x_, shift, scale
    else:
        return x_


def pseudoinverse(M, return_condition_number=False):
    '''Find the pseudoinverse of the matrix M using singular value
       decomposition.
    INPUT
        M - array_like
            An (m x n) matrix whose pseudoinverse we want to find.
        return_condition_number - bool [default: False]
            If True, the function will also return the condition number
            associated with the inversion.
    OUTPUT
        pinv - array_like
            The Moore-Penrose pseudoinverse of the matrix M.
        cond - float
            The condition number of the inversion (i.e., the ratio of the
            largest to smallest singular values).
    '''
    # Compute the singular value decomposition.
    U, s, Vt = np.linalg.svd(M, full_matrices=False)
    V = Vt.T

    # Construct the pseudoinverse and compute the condition number.
    M_pinv = V.dot(np.diag(1/s)).dot(U.T)
    condition_number = s.max()/s.min() 

    # If requested, return condition number; otherwise, don't.
    if return_condition_number:
        return M_pinv, condition_number
    else:
        return M_pinv


def least_squares(basis, y=None, coefs=None):
    '''Find the least square fit to one-dimensional data using the specified 
       basis. OR, if coefficients are specified rather than data, y, the
       coefficients are used to generate the fit; in this case, the number of
       coefficients must equal basis.nb_bases.
    INPUT
        basis - basis object
            A mathtools basis object (see mathtools.legendre, e.g.)
        y - array_like
            One dimensional array of data to fit.
        coefs - array_like
            Coefficients used to project basis generate fit.
    OUTPUT
       fit - Struct
            A Struct containing the following fields:
                - x:        The domain on which fit is defined.
                - y:        The best fit to the data.
                - dy:       The derivative of the best fit.
                - d2y:      The second derivative of the best fit.
                - coefs:    The coefficients used for the fit.
    '''

    # Do we have data, or must we do the fit?
    if (y is not None):
        # Augment the data (for regularization).
        augmented_y = basis.augment(y)

        # Perform the least squares fit using the precomputed pseudoinverse!
        coefs = basis.inverse.dot(augmented_y)

    if (coefs is None):
        raise ValueError('Cannot project onto basis! No data or coefficients!')

    # Generate a fit Struct object to hold the results.
    fit                         = Struct()
    fit.coefs                   = coefs
    fit.x                       = basis.x
    fit.y                       = np.zeros_like(basis.x)
    fit.dy                      = np.zeros_like(basis.x)
    fit.d2y                     = np.zeros_like(basis.x)
    fit.y[basis.valid_idx]      = basis.B.dot(coefs)
    fit.dy[basis.valid_idx]     = basis.dB.dot(coefs)
    fit.d2y[basis.valid_idx]    = basis.d2B.dot(coefs)

    return fit



class Vessel(object):
    '''Create a container object that holds properties. Can be easily saved &
       loaded.  
    USAGE
        Create a Vessel instance:
            >>> data = Vessel('storage.dat')
        Assign properties to the object:
            >>> data.variable = [1,2,3,'string data']
            >>> data.my_array = np.arange(0,100,0.1)
            >>> data.my_dict = {'name': 'Albert Einstein'}
        Save the object! If no filename is specified, it will use the initially
        supplied filename, if present; otherwise, an error is thrown.
            >>> data.save()
        When we want to load the data, simply create another instance. A
        filename may be passed during object creation. If the filename
        corresponds to an existing file, the file will automatically be loaded:
            >>> other = Vessel('storage.dat')
        Otherwise, a file may be loaded explicitly, at some later point:
            >>> other.load('other_file.dat')
            >>> other.variable # ==> [1,2,3,'string data']
            >>> other.my_dict  # ==> {'name': 'Albert Einstein'}
            >>> other.my_array # ==> array([ 0. ,  0.1,  0.2,  0.3, ... 9.9])
            >>> other.keys     # ==> ['my_array', 'my_dict', 'variable',
            >>>                #      'current_filename']
        When the .save() method is later called, the current filename will be
        used, unless another filename is explicitly specified as a parameter to
        the save command:
            >>> other.save('new_file.dat') # ==> Saved to a new file!
    TIPS
        To determine the properties attached to an object instance, examine the
        .keys property. This will list the names of all attributes attached to
        the instance.
    INGESTING DICTIONARIES
        The Vessel object also allows for the ingestion of large dictionaries.
        This is useful for saving all variables in the local namespace. As an
        example:
            >>> ignore_vars = locals() 
            >>> x = 42; y = np.sin(pi/3); z = np.arange(0,5,0.1)
            >>> v = Vessel('kitchen_sink.data)
            >>> v.ingest(locals(), ignore_vars)
            >>> v.save()

        We have now grabbed all variables from the local scope and saved them
        to disk. We can reconstitute this scope at a later time as follows:

            >>> w = Vessel('kitchen_sink.data') # loads data if the file exists
            >>> for key in w.keys:
            >>>     exec('%s=w.%s') % (key,key)

        The previously saved local scope will now be reconstituted.
    '''

    def __init__(self, filename=None):
        self._filename = filename
        if self._filename:
            # If filename specified, and file exists, load it.
            if len(glob.glob(filename)) > 0:
                self.load()

    def _set_filename(self, filename):
        '''Set the object's filename. If filename does not exist, throw an
        error.'''
        if filename:
            self._filename = filename
        if not self._filename:
            raise ValueError('No filename specified.')

    def ingest(self, var_dict, ignore_variable_names=None):
        '''Ingest a dictionary of variables (such as locals(), e.g.). Only
        variables in the supplied (or default) white list will be retained.
        Variables are added as attributes to the object.  '''
        if ignore_variable_names:
            self.ignore_variable_names = ignore_variable_names
        else:
            self.ignore_variable_names = []
        for key in var_dict.keys():
            if key not in self.ignore_variable_names:
                self.__dict__[key] = var_dict[key]

    @property
    def keys(self):
        keys = list(self.__dict__.keys())
        keys.remove('_filename') # don't show internal filename.
        keys.sort()
        keys.append('current_filename')
        return keys

    @property
    def current_filename(self):
        return self._filename

    def save(self, filename=None):
        '''Save the data into a file with the specified name.'''
        self._set_filename(filename)
        f = open(self._filename, 'wb')
        pickle.dump(self.__dict__, f, protocol=pickle.HIGHEST_PROTOCOL)
        f.close()

    def load(self, filename=None):
        '''Load object from specified file.'''
        self._set_filename(filename)
        f = open(self._filename, 'rb')
        loaded_object = pickle.load(f)
        f.close()
        # Unpack the object and add variables as properties to this object.
        for key, val in loaded_object.items():
            self.__dict__[key] = val


def mahal(x, mu=None, S=None, return_stats=False):
    '''Find the Mahalanobis distance between row vectors in x and mean mu, with
       covariance S. If mu and S are not provided, the mean and covariance of
       x are used instead.
    INPUTS
        x - array_like
            A matrix of row vectors.
        mu - array_like
            Mean vector. mu.shape[1] must equal x.shape[1].
        S - array_like
            The covariance matrix. S.shape[0] == S.shape[1] == x.shape[1]
        return_stats - boolean
            If True, returns the mean and covariance used in the calcuations.
    OUTPUTS
        mahal_dist - array_like
            An array of the Mahalanobis distances associated with the vectors
            in x.
        mu - array_like
            The mean vector used in Mahalanobis calculations.  
        S - array_like
            The covariance matrix used.
        '''
    x = np.array(x)
    if (mu is None):
        mu = x.mean(0)
        S = np.cov(x.T)

    mahal_dist = np.zeros(x.shape[0])
    inv_S = np.linalg.pinv(S)
    for i, x_i in enumerate(x):
        mahal_dist[i] = (x_i - mu).T.dot(inv_S).dot((x_i - mu))

    if return_stats:
        return (mahal_dist, mu, S)
    else:
        return mahal_dist
