# Math Tools

The ```mathtools```  package provides algorithms for processing and
manipulating data. For more details on the available tools, please [see
below](#tools).

## Installation

Math tools is currently installed by cloning the git repository to your local
machine and adding the ```mathtools``` path to your Python path. To
clone this repository, simply:

```unix
$~> git clone https://github.com/michiganaerospace/mathtools.git
```

Enter your git username and password, and the repository will be cloned into
your local directory, say ```/Users/mjl/develop/local_python/mathtools```. In
your ```.bash_rc``` file (or similar), add this path to your ```PYTHONPATH```:

```unix
export PYTHONPATH="${PYTHONPATH}:/Users/mjl/develop/local-python/mathtools"
```

Now you can import the ```mathtools``` package anywhere on your machine:

```python
from mathtools.fit import Fit

# Create a fit object on points x:
f = Fit(x, nb_bases=15)
```

## Tests

All code in the ```mathtools``` is fully tested. The tests are available in the
[tests](https://github.com/michiganaerospace/mathtools/tree/master/tests)
directory, and can be run using the nose testing framework. For a good
introduction to nose, check out this [resource](http://goo.gl/hfzYPz). 

As a general rule, before any code is checked in to the repository, one should
verify that all unit tests are passing.

## Tools 

- ```fit```
    - [```fit.Fit```](#fit) — The ```Fit``` class provides algorithms for
      regularized least squares using different bases.
    - [```fit.best_fit```](#best_fit) — Perform a least square fit
      to data, given a basis object.
- ```fourier```
    - [```fourier.fourier_basis```](#fourier) — Function for generating a
      Fourier series basis.
    - [```fourier.d_fourier_basis```](#d_fourier) — Function for generating the
      derivatives of Fourier series basis.
    - [```fourier.d2_fourier_basis```](#d2_fourier) — Function for generating
      the second derivatives of a Fourier series basis.
    - [```fourier.create_fourier_basis```](#create_fourier) — Generate a
      Fourier basis object.
- ```legendre```
    - [```legendre.legendre_basis```](#legendre) — Function for generating a
      Legendre polynomial basis.
    - [```legendre.d_legendre_basis```](#d_legendre) — Generates the first
      derivative of a Legendre polynomial basis.
    - [```legendre.d2_legendre_basis```](#d2_legendre) — Generates the second
      derivative of a Legendre polynomial basis.
    - [```legendre.create_legendre_basis```](#create_legendre) — Generate a
      Legendre basis object.
- ```splines```
    - [```splines.cubic_spline_basis_knot_interval```](#spline_knot_basis)
      — Generates a cubic spline basis on the unit interval.
    - [```splines.d_cubic_spline_basis_knot_interval```](#spline_knot_d_basis)
      — Generates the derivative of the cubic spline basis between knots. 
    - [```splines.d2_cubic_spline_basis_knot_interval```](#spline_knot_d2_basis)
      — Generates the derivative of the cubic spline basis between knots.
    - [```splines.cubic_spline_basis```](#spline_basis)
      — Generates the cubic spline basis on the domain.
    - [```splines.d_cubic_spline_basis```](#d_spline_basis)
      — Generates the derivative of the cubic spline basis on the domain.
    - [```splines.d2_cubic_spline_basis```](#d2_spline_basis)
      — Generates the derivative of the cubic spline basis on the domain.
    - [```splines.create_spline_basis```](#create_spline_basis)
      — Creates a cubic spline basis object. 
- ```utils```
    - [```utils.map_to_interval```](#map_interval) — Map an array into
      specified interval.
    - [```utils.pseudoinverse```](#pseudoinverse) — Find the pseudoinverse of a
      matrix M using singular value decomposition.
    - [```utils.Vessel```](#vessel) — Provides a convenient container for
      storing, saving, and loading data.

<a name='fit'></a>
### ```fit.Fit```

The ```Fit``` class is a machine that allows you to easily fit data using
regularized least squares. It provides a convenient interface to the core basis
generation and fitting routines. To see how to use this class, jump right to 
the [usage](#fit_usage) section. Additionally, example code can be found in the
[scripts](https://goo.gl/tQl9E8) directory.

> **METHODS**
   
> **```Fit(x=None, nb_bases=0, basis_type='legendre', freq=1.0, reg_coefs=[0,0,0])```**

> Creates a new ```Fit``` object.   

> > ARGUMENTS
> > - **```x — array_like```**: An array of abscissa values — an
> >   ```nb_samples``` length vector of 'x values'.
> > - **```nb_bases — int```**: The number of basis vectors to use when fitting
> >   the data.  In the case of a cubic spline basis, this corresponds to the
> >   number of knots used.
> > - **```basis_type — str```**: The type of basis to use for the fitting. May
> >   have values of ```legendre```, ```fourier```, or ```cubic-spline```.
> > - **```freq — float```**: The fundamental frequency of a Fourier series
> >   basis. If a basis other than Fourier is selected, this parameter is
> >   ignored. 
> > - **```reg_coefs — array_like```**: A list or array of three regularization
> >   coefficients for penalizing the magnitude of the fit and its first and
> >   second derivatives, respectively. The default value is ```reg_coefs=[0.0,
> >   0.0, 0.0]```.   
>
> > OUTPUT    
> > - None 
> 
> **```Fit.fit(y)```**
>
> The method fits the data, ```y```, using the current configurations of the
> ```Fit``` object.

> > ARGUMENTS   
> > - **```y — array_like```**: Vector of data samples that we want to fit. The
> >   data in the vector ```y``` must correspond to the domain ```x``` for
> >   which the basis was constructed.
> >
> > OUTPUT  
> >   - **```results — object```**: An object containing fit results. It has
> >     the following properties:
> >       - ```x — array_like```: The domain associated with the fit.
> >       - ```y — array_like```: The fit sampled on the domain. 
> >       - ```dy — array_like```: The derivative of the fit.
> >       - ```d2y — array_like```: The second derivative of the fit.
> >       - ```coefs — array_like```: The coefficients of the fit.
> 
> **```Fit.resample(x)```**   
>
> Resamples the current fit to the specified domain ```x``` using the existing
> coefficients. If any part of ```x``` lies outside the original domain of the
> fit, the values for the fit in that region are set to zero. The ```fit()```
> must have been run before the fit can be resampled.
>
> > ARGUMENTS
> >   - **```x — array_like```**: A vector of abscissa values — an
> >     ```nb_samples``` length vector of x-values.
> > 
> > OUTPUT
> >  - **```results — object```**: An object containing fit results. It has the
> >    following properties:
> >      - ```x — array_like```: The domain associated with the fit.
> >      - ```y — array_like```: The fit sampled on the domain.
> >      - ```dy — array_like```: The derivative of the fit.
> >      - ```d2y — array_like```: The second derivative of the fit.
> >      - ```coefs — array_like```: The coefficients of the fit.
>
> **```Fit.config(x=None, nb_bases=0, basis_type=None, freq=1.0, reg_coefs=None)```**
> 
> Recomputes the basis given a change in the underlying parameters. When the
> basis is updated, all existing fit coefficients are discarded. 
> 
> > ARGUMENTS    
> > - Defined as [above](#fit).
> >
> > OUTPUT   
> >   - None.
>
> **PROPERTIES** <a name='fit-properties'></a>
> 
> When the ```Fit``` class is instantiated, it makes several properties
> available to the user.
>
>   - **```basis — object```**: The current basis object. For details on what
>     basis objects are and how to use them, check out the discussion of [basis
>     objects](#basis_object).
>   - **```basis_type — str```**: The type of basis to use for the fitting. May
>     have values of ```legendre```, ```fourier```, or ```cubic-spline```.
>   - **```coefs — array_like```**: The current fit coefficients. If there is
>     no active fit, the coefficients are set to  ``None``.
>   - **```nb_bases — int```**: The number of basis vectors currently used in
>     the basis.
>   - **```reg_coefs — array_like```**: A list or array of three regularization
>     coefficients for penalizing the magnitude of the fit and its first and
>     second derivatives, respectively. 
>   - **```x — array_like```**: The domain on which the basis was built.
>
> <a name='fit_usage'></a>
> **USAGE** 
> 
> To get a sense of what the ```Fit``` class can do, let's try to fit some
> noisy data. We'll generate a sine wave and add a little noise to it, as
> illustrated in the following code snippet.
>
> ```python 
> # Create some noisy data.  
> t = np.linspace(0,15*np.pi, 300) 
> y = np.sin(2*np.pi/5*t) + 0.2 * np.random.randn(len(t)) 
> ```
> 
> The noisy sine data is shown in the following figure.
> 
> <img src='docs/images/noisy_sine.png' width='500'/>
> 
> We'd like to fit this data with a smooth curve. We can use the ```Fit```
> class to do this quite easily. We first create an instance of the ```Fit```
> class, intializing it with the domain on which we'd like to fit the data, and
> the number of basis vectors we want to use, 
> 
> ```python 
> # Create a fit object with 15 basis vectors.  
> f = Fit(t, 15) 
> ```
> 
> This will create an instance of a ```Fit``` object, generate a Legendre
> polynomial basis for the provided abscissa (here, the vector ```t```) with 15
> basis vectors. The Legendre polynomial basis is the default; others may be
> specified. To fit the data, we use the object's ```fit``` method,
> 
> ```python 
> # Fit the noisy data.  
> r = f.fit(y) 
> ```
> 
> The fit method returns a structure, here ```r```, which contains the fit
> coefficients, the fit itself, and its first two derivatives. For convenience,
> the original domain is also included in the structure. To see how well we
> did, let's plot the fit on top of the original, noisy data.
> 
> ```python 
> plot(r.x, r.y, '-.', linewidth=2) 
> ```
> 
> <img src='docs/images/noisy_sine_fit.png' width='500'/>
> 
> Note that although the blue curve smoothly fits the red data points, we are
> not explicitly penalizing the derivatives here (the regularization
> coefficients are, by default, zero). We are instead *effectively*
> regularizing by using a small number of Legendre polynomial basis vectors —
> in this case, fifteen basis vectors. For more on this topic, check out the
> technical discussion at ...
> 
> But what if we wanted to sample this function at a different set of points?
> As it stands, our smooth curve is sampled only at the values of ```x```
> corresponding to the original data. If we wish to sample our smooth curve on
> a different set of data points, we can use the ```Fit``` object's
> ```resample``` method. For example, suppose we wish to sample on a larger
> number of points, and on a different domain.
> 
> ```python 
> t_new = np.linspace(-np.pi, 6*np.pi, 500) 
> ```
> 
> Note that we are now sampling at 500 points, extending from -pi to 6pi. Using
> the previously computed ```Fit``` object, we compute a new results object,
> 
> ```python 
> rs = f.resample(t_new) 
> ```
> 
> The new fit now samples the same curve (i.e., using the same fit
> coefficients), but at different points. In this example, we have deliberately
> resampled the curve on a domain that extends beyond the support of the
> original data. Of course, we cannot expect to fit data in this region (i.e.
> the ```Fit``` object does not extrapolate). The convention here is to set the
> fit to zero in those regions that do not intersect the original domain. We
> illustrate this by plotting the resampled data:
> 
> ```python 
> plot(rs.x, rs.y, '-.', linewidth=2) 
> ```
> 
> <img src='docs/images/noisy_sine_new_domain.png' width='500'/>
> 
> Where the resampled domain intersects the support of the original data, we
> reproduce the fit. However, once we venture beyond the support of that data,
> the fit returns zero.
> 
> As mentioned above, the ```results``` object returned from the ```fit```
> method also contains the derivatives of the fit, which may be useful. To take
> a look at the first derivative, for example, we can plot ```r.dy```:
> 
> ```python 
> plot(r.x, r.dy, '-.', linewidth=2) 
> ```
> 
> <img src='docs/images/noisy_sine_fit_deriv.png' width='500'/>
> 
> We may want to alter our fit in some way — perhaps by choosing a different
> basis, or selecting a different number of basis vectors. The ```Fit```
> class's ```config``` method allows us to easily change parameters and
> recompute the basis. For example, suppose we wished to use only ten basis
> vectors, rather than fifteen, as above. We simply call ```config()``` with
> the parameters we wish to change:
> 
> ```python 
> f.nb_bases # ==> 15 
> f.config(nb_bases=10) 
> f.nb_bases # ==> 10 
> ```
> 
> Other parameters can be changed similarly. When ```config``` is called, the
> basis is recomputed and all fit coefficients are discarded.

<a name='best_fit'></a>
### ```fit.best_fit(basis, y=None, coefs=None)```

Find the least square fit to one-dimensional data using the specified basis. If
coefficients are specified rather than data, the coefficients are used to
generate the fit. The function returns a structure, described in detail below,
which provides the fit, its derivatives, and the coefficients.

This function powers the ```Fit``` class. 

> ARGUMENTS
>   - **```basis```**: A basis object. See the detailed discussion
>     [here](#basis_object).
>   - **```y — array_like```**: An ```nb_samples``` length data array that we
>     wish to fit. Providing ```y``` is optional; if it is not provided, fit
>     coefficients must be provided.
>   - **```coefs — array_like```**: Coefficients defining the fit. This is
>     useful for resampling a curve using a different basis. The number of
>     coefficients must be compatible with the basis object. ```coefs``` is
>     optional; but if it is not provided, data must be provided.
> 
> OUTPUT
>   - **```fit```**: A fit object with the following properties: 
>       - ```x — array_like```: The domain on which the fit is defined.
>       - ```y — array_like```: The best fit to the data.
>       - ```dy — array_like```: The derivative of the best fit to the data.
>       - ```d2y — array_like```: The second derivative of the best fit to the
>         data.
>       - ```coefs — array_like```: The coefficients used for the fit.
>
> USAGE
> ```python
> import numpy as np
> from mathtools.legendre import create_legendre_basis
> from mathtools.fit import best_fit
> 
> # Generate some noisy data to fit.
> x = np.linspace(0, 3*np.pi, 200)
> y = np.cos(x/4) + np.random.randn(len(x))
> 
> # Create a basis for fitting this noisy data.
> basis = create_legendre_basis(x, 15)
> 
> # Fit the noisy data!
> fit = best_fit(basis, y)
> 
> # Check out the coefficients.
>  fit.coefs # ==>
> # array([ 0.34114583, -0.95155628, -0.47158956,  0.03764087,  0.14931675,
> #        -0.27857014, -0.26873534,  0.28205916,  0.64860036, -0.16132168,
> #        -0.11036244,  0.13654455,  0.03552914,  0.14698555,  0.33255431])
> ```

<a name='fourier'></a>
### ```fourier.fourier_basis(x, nb_bases, freq=1.0)```

Computes a Fourier series basis on the specified domain.

> ARGUMENTS    
>   - **```x — array_like```**: The domain over which we are defining the
>     basis. An ```nb_samples``` length vector.
>   - **```nb_bases — int```**: The number of basis vectors to generate.
>   - **```freq — float```**: The fundamental frequency of the series.
>
> OUTPUT   
>   - **```B — array_like```**: An ```nb_samples x nb_bases``` array. The
>     columns of the array are the Fourier series basis vectors.


<a name='d_fourier'></a>
### ```fourier.d_fourier_basis(x, nb_bases, freq=1.0)```

Computes a basis corresponding to the first derivative of a Fourier series on
the specified domain.

> ARGUMENTS    
>   - **```x — array_like```**: The domain over which we are defining the
>     basis. An ```nb_samples``` length vector.
>   - **```nb_bases — int```**: The number of basis vectors to generate.
>   - **```freq — float```**: The fundamental frequency of the series.
>
> OUTPUT   
>   - **```dB — array_like```**: An ```nb_samples x nb_bases``` array. The
>     columns of the array are the derivatives of the Fourier series basis
>     vectors.


<a name='d2_fourier'></a>
### ```fourier.d2_fourier_basis(x, nb_bases, freq=1.0)```

Computes a basis corresponding to the second derivative of a Fourier series on
the specified domain.

> ARGUMENTS    
>   - **```x — array_like```**: The domain over which we are defining the
>     basis. An ```nb_samples``` length vector.
>   - **```nb_bases — int```**: The number of basis vectors to generate.
>   - **```freq — float```**: The fundamental frequency of the series.
>
> OUTPUT   
>   - **```d2B — array_like```**: An ```nb_samples x nb_bases``` array. The
>     columns of the array are the second derivatives of the Fourier series
>     basis vectors.



<a name='create_fourier'></a>
### ```fourier.create_fourier_basis(x, nb_bases, freq=1.0, reg_coefs=[0,0,0], x_ref=None)```

This function creates a Fourier series basis object. The basis object contains
everything required to perform a regularized least squares fit of data on the
specified domain. Basis objects are typically used to fit data through the use
of the [```best_fit```](#best_fit) helper function, or as part of the
[```Fit```](#fit) class.

> ARGUMENTS    
>   - **```x — array_like```**: The domain over which we are defining the
>     basis. An ```nb_samples``` length vector.
>   - **```nb_bases — int```**: The number of basis vectors to generate.
>   - **```freq — float```**: The fundamental frequency of the series.
>   - **```reg_coefs — array_like```**: A list or array of three regularization
>     coefficients for penalizing the magnitude of the fit and its first and
>     second derivatives, respectively. 
>   - **```x_ref — array_like```**: An optional reference domain. This is
>     useful for resampling data.  It ensures that data is mapped to the
>     interval [0, 1] in a consistent manner, and allows us to avoid
>     attempting to fit data outside of the original domain.
>
> OUTPUT <a name='basis_object'></a>   
>   - **```basis — object```**: A basis object. Basis objects consolidate all
>     the information required to fit data. The basis object contains the
>     following properties and methods:
>       - **```augment(y)```**: A method that takes in an ```nb_samples```
>         length data vector, ```y```, and returns a properly augmented data
>         vector.  The vector is concatenated with the proper number of zeros
>         so that regularized least squares just works.
>       - **```x```**: The domain over which the basis is defined.
>       - **```valid_idx```**: The indices of valid entries in the domain. If a
>         reference domain, ```x_ref``` is specified, regions of the domain
>         ```x``` that fall outside of the reference domain are considered
>         invalid, and are excluded when the fit is computed.
>       - **```B — array_like```**: An ```nb_samples x nb_bases``` array of
>         Legendre polynomial basis (column) vectors.
>       - **```B_ — array_like```** The so-called *brick*. The brick is a
>         vertical concatenation of the ```B```, ```I```, ```dB```, ```d2B```
>         matrices.  The brick matrix allows us to force the solution and its
>         derivative to be close to zero, as dictated by the regularization
>         coefficients.  In fact, only those components with nonzero
>         regularization coefficients are included in the brick, in order to
>         minimize computational overhead during computation of the SVD. The
>         matrix ```I``` is an ```nb_bases x nb_bases``` sized identity
>         matrix. It serves to penalize the L2 norm of the fit coefficients.  
>       - **```condition_number — float```**: The condition number associated
>         with the pseudoinversion of the matrix ```B_```.
>       - **```inverse — array_like```**: The pseudoinverse of the ```B_```
>         matrix. May be used to compute the fit coefficients to a data vector.
>         For example, to find the fit coefficients to a vector ```y```, we
>         compute ```basis.inverse.dot(basis.augment(y))```, where we have also
>         used the ```basis.augment()``` method to ensure the data vector has
>         the appropriate dimensions.
>       - **```nb_bases — int```**: Number of basis vectors used.
>       - **```reg_coefs```**: The regularization coefficients that were used
>         to construct the basis.
>       - **```dB — array_like```**: Derivative of basis vectors in ```B```.
>         An ```nb_samples x nb_bases``` sized array.
>       - **```d2B — array_like```**: Second derivative of basis vectors in
>         ```B```. An ```nb_samples x nb_bases``` sized array.

<a name="legendre"></a>
### ```legendre.legendre_basis(x, nb_bases)```

Computes the Legendre polynomial basis on a specified domain.

> ARGUMENTS    
>   - **```x — array_like```**: The domain over which we are defining the
>     basis. An ```nb_samples``` length vector.
>   - **```nb_bases — int```**: The number of basis vectors to generate.
>
> OUTPUT   
>   - **```B — array_like```**: An ```nb_samples x nb_bases``` array. The
>     columns of the array are the Legendre polynomial basis vectors.

<a name="d_legendre"></a>
### ```legendre.d_legendre_basis(x, nb_bases)```

Computes the derivative of the Legendre polynomial basis on a specified domain.

> ARGUMENTS    
>   - **```x — array_like```**: The domain over which we are defining the
>     basis. An ```nb_samples``` length vector.
>   - **```nb_bases — int```**: The number of basis vectors to generate.
>
> OUTPUT   
>   - **```dB — array_like```**: An ```nb_samples x nb_bases``` array. The
>     columns of the array are the first derivatives of the Legendre polynomial
>     basis vectors.


<a name="d2_legendre"></a>
### ```legendre.d2_legendre_basis(x, nb_bases)```

Computes the second derivative of the Legendre polynomial basis on a specified
domain.

> ARGUMENTS    
>   - **```x — array_like```**: The domain over which we are defining the
>     basis. An ```nb_samples``` length vector.
>   - **```nb_bases — int```**: The number of basis vectors to generate.
>
> OUTPUT   
> - **```d2B — array_like```**: An ```nb_samples x nb_bases``` array. The
>   columns of the array are the second derivativeis of the Legendre polynomial
>   basis vectors.


<a name='create_legendre'></a>
### ```legendre.create_legendre_basis(x, nb_bases, reg_coefs=[0,0,0], x_ref=None)```

This function creates a Legendre polynomial basis object. The basis object 
contains everything required to perform a regularized least squares fit of
data on the specified domain. Basis objects are typically used to fit data
through the use of the [```best_fit```](#best_fit) helper function,
or as part of the [```Fit```](#fit) class.

> ARGUMENTS    
>   - **```x — array_like```**: The domain over which we are defining the
>     basis. An ```nb_samples``` length vector.
>   - **```nb_bases — int```**: The number of basis vectors to generate.
>   - **```reg_coefs — array_like```**: A list or array of three regularization
>     coefficients for penalizing the magnitude of the fit and its first and
>     second derivatives, respectively. 
>   - **```x_ref — array_like```**: An optional reference domain. This is
>     useful for resampling data.  It ensures that data is mapped to the
>     interval [-1, 1] in a consistent manner, and allows us to avoid
>     attempting to fit data outside of the original domain.
>
> OUTPUT <a name='basis_object'></a>   
>   - **```basis — object```**: A basis object. Basis objects consolidate all
>     the information required to fit data. The basis object contains the
>     following properties and methods:
>       - **```augment(y)```**: A method that takes in an ```nb_samples```
>         length data vector, ```y```, and returns a properly augmented data
>         vector.  The vector is concatenated with the proper number of zeros
>         so that regularized least squares just works.
>       - **```x```**: The domain over which the basis is defined.
>       - **```valid_idx```**: The indices of valid entries in the domain. If a
>         reference domain, ```x_ref``` is specified, regions of the domain
>         ```x``` that fall outside of the reference domain are considered
>         invalid, and are excluded when the fit is computed.
>       - **```B — array_like```**: An ```nb_samples x nb_bases``` array of
>         Legendre polynomial basis (column) vectors.
>       - **```B_ — array_like```** The so-called *brick*. The brick is a
>         vertical concatenation of the ```B```, ```I```, ```dB```, ```d2B```
>         matrices.  The brick matrix allows us to force the solution and its
>         derivative to be close to zero, as dictated by the regularization
>         coefficients.  In fact, only those components with nonzero
>         regularization coefficients are included in the brick, in order to
>         minimize computational overhead during computation of the SVD. The
>         matrix ```I``` is an ```nb_bases x nb_bases``` sized identity
>         matrix. It serves to penalize the L2 norm of the fit coefficients.  
>       - **```condition_number — float```**: The condition number associated
>         with the pseudoinversion of the matrix ```B_```.
>       - **```inverse — array_like```**: The pseudoinverse of the ```B_```
>         matrix. May be used to compute the fit coefficients to a data vector.
>         For example, to find the fit coefficients to a vector ```y```, we
>         compute ```basis.inverse.dot(basis.augment(y))```, where we have also
>         used the ```basis.augment()``` method to ensure the data vector has
>         the appropriate dimensions.
>       - **```nb_bases — int```**: Number of basis vectors used.
>       - **```reg_coefs```**: The regularization coefficients that were used
>         to construct the basis.
>       - **```dB — array_like```**: Derivative of basis vectors in ```B```.
>         An ```nb_samples x nb_bases``` sized array.
>       - **```d2B — array_like```**: Second derivative of basis vectors in
>         ```B```. An ```nb_samples x nb_bases``` sized array.
 

<a name='spline_knot_basis'></a>
### ```splines.cubic_spline_basis_knot_interval(x)```

Compute the four basis vectors over a knot interval for x.

> ARGUMENTS
>   - **```x — array_like```**: A one dimensional array of data; the domain on
>     which we wish to fit the data.
> 
> OUTPUT
>   - **```B — array_like```**: The four basis vectors defining the cubic
>     spline on the unit interval — an (nb_samples x 4) size array. 


<a name='spline_knot_d_basis'></a>
### ```splines.d_cubic_spline_basis_knot_interval(x)```

Compute the four derivative basis vectors over the knot interval for x.

> ARGUMENTS
>   - **```x — array_like```**: A one dimensional array of data; the domain on
>     which we wish to fit the data.
> 
> OUTPUT
>   - **```dB — array_like```**: The four basis vectors defining the derivative
>     of the cubic spline on the unit interval — an (nb_samples x 4) size
>     array.

 
<a name='spline_knot_d2_basis'></a>
### ```splines.d2_cubic_spline_basis_knot_interval(x)```

Compute the four second derivative basis vectors over the knot interval for x.

> ARGUMENTS
>   - **```x — array_like```**: A one dimensional array of data; the domain on
>     which we wish to fit the data.
> 
> OUTPUT
>   - **```dB — array_like```**: The four basis vectors defining the derivative
>     of the cubic spline on the unit interval — an (nb_samples x 4) size
>     array. 


<a name='spline_basis'></a>
### ```splines.cubic_spline_basis(x, nb_knots)```

Compute the cubic spline basis for domain ```x``` using nb_knots uniformly
distributed knots.

> ARGUMENTS
>   - **```x — array_like```**: A one dimensional array of data; the domain on
>     which we wish to fit the data.
>   - **```nb_knots — array_like```**: The number of uniformly distributed
>     knots to use in constructing the basis.
> 
> OUTPUT
>   - **```B — array_like```**: The ```nb_samples x nb_knots``` array
>     containing (column) cubic spline basis vectors.


<a name='d_spline_basis'></a>
### ```splines.d_cubic_spline_basis(x, nb_knots)```

Compute the derivative of the cubic spline basis for domain ```x``` using
nb_knots uniformly distributed knots.

> ARGUMENTS
>   - **```x — array_like```**: A one dimensional array of data; the domain on
>     which we wish to fit the data.
>   - **```nb_knots — array_like```**: The number of uniformly distributed
>     knots to use in constructing the basis.
> 
> OUTPUT
>   - **```dB — array_like```**: The ```nb_samples x nb_knots``` array
>     containing (column) derivative cubic spline basis vectors.


<a name='d2_spline_basis'></a>
### ```splines.d2_cubic_spline_basis(x, nb_knots)```

Compute the second derivative of the cubic spline basis for domain ```x```
using nb_knots uniformly distributed knots.

> ARGUMENTS
>   - **```x — array_like```**: A one dimensional array of data; the domain on
>     which we wish to fit the data.
>   - **```nb_knots — array_like```**: The number of uniformly distributed
>     knots to use in constructing the basis.
> 
> OUTPUT
>   - **```d2B — array_like```**: The ```nb_samples x nb_knots``` array
>     containing (column) second derivative cubic spline basis vectors.

<a name='create_spline_basis'></a>
### ```splines.create_spline_basis(x, nb_bases, reg_coefs=[0,0,0], x_ref=None)```

This function creates a cubic spline basis object. The basis object contains
everything required to perform a regularized least squares fit of data on the
specified domain. Basis objects are typically used to fit data through the use
of the [```best_fit```](#best_fit) helper function, or as part of the
[```Fit```](#fit) class.

> ARGUMENTS    
>   - **```x — array_like```**: The domain over which we are defining the
>     basis. An ```nb_samples``` length vector.
>   - **```nb_bases — int```**: The number of basis vectors to generate.
>   - **```reg_coefs — array_like```**: A list or array of three regularization
>     coefficients for penalizing the magnitude of the fit and its first and
>     second derivatives, respectively. 
>   - **```x_ref — array_like```**: An optional reference domain. This is
>     useful for resampling data.  It ensures that data is mapped to the
>     interval [0, 1] in a consistent manner, and allows us to avoid
>     attempting to fit data outside of the original domain.
>
> OUTPUT <a name='basis_object'></a>   
>   - **```basis — object```**: A basis object. Basis objects consolidate all
>     the information required to fit data. The basis object contains the
>     following properties and methods:
>       - **```augment(y)```**: A method that takes in an ```nb_samples```
>         length data vector, ```y```, and returns a properly augmented data
>         vector.  The vector is concatenated with the proper number of zeros
>         so that regularized least squares just works.
>       - **```x```**: The domain over which the basis is defined.
>       - **```valid_idx```**: The indices of valid entries in the domain. If a
>         reference domain, ```x_ref``` is specified, regions of the domain
>         ```x``` that fall outside of the reference domain are considered
>         invalid, and are excluded when the fit is computed.
>       - **```B — array_like```**: An ```nb_samples x nb_bases``` array of
>         Legendre polynomial basis (column) vectors.
>       - **```B_ — array_like```** The so-called *brick*. The brick is a
>         vertical concatenation of the ```B```, ```I```, ```dB```, ```d2B```
>         matrices.  The brick matrix allows us to force the solution and its
>         derivative to be close to zero, as dictated by the regularization
>         coefficients.  In fact, only those components with nonzero
>         regularization coefficients are included in the brick, in order to
>         minimize computational overhead during computation of the SVD. The
>         matrix ```I``` is an ```nb_bases x nb_bases``` sized identity
>         matrix. It serves to penalize the L2 norm of the fit coefficients.  
>       - **```condition_number — float```**: The condition number associated
>         with the pseudoinversion of the matrix ```B_```.
>       - **```inverse — array_like```**: The pseudoinverse of the ```B_```
>         matrix. May be used to compute the fit coefficients to a data vector.
>         For example, to find the fit coefficients to a vector ```y```, we
>         compute ```basis.inverse.dot(basis.augment(y))```, where we have also
>         used the ```basis.augment()``` method to ensure the data vector has
>         the appropriate dimensions.
>       - **```nb_bases — int```**: Number of basis vectors used.
>       - **```reg_coefs```**: The regularization coefficients that were used
>         to construct the basis.
>       - **```dB — array_like```**: Derivative of basis vectors in ```B```.
>         An ```nb_samples x nb_bases``` sized array.
>       - **```d2B — array_like```**: Second derivative of basis vectors in
>         ```B```. An ```nb_samples x nb_bases``` sized array.

<a name='map_interval'></a>
### ```utils.map_to_interval(x, interval, return_all=False)```

Shift and scale array so that its elements live in the specified interval. The
function finds shift and scale factors,

```python
x_ = scale * (x - shift)
```

such that the elements of ```x_``` are in the interval specified.  If
requested, the function will return the scale and shift factors that perform
this mapping.

> ARGUMENTS
>   - **```x — array_like```**: A one dimensional array of data. 
>   - **```interval — array_like```**: A two element list, array, or tuple
>     defining the minimum and maximum boundaries of the interval into which we
>     wish to shift and scale the data in the array ```x```.
>   - **```return_all — Boolean```**: Boolean flag (defaults to false), which
>     determines whether to return the shift and scale factors for later use.
> 
> OUTPUT
>   - **```x_ — array_like```**: The original data array, scaled and shifted so
>     that its elements live in the specified interval.
>   - **```shift — float```**: The shift factor. Only returned if
>     ```return_all``` is set to ```True```.  
>   - **```scale```**: The scale factor. Only returned if ```return_all``` is
>     set to ```True```.  


<a name='pseudoinverse'></a>
### ```utils.pseudoinverse(M, return_condition_number=False)```

Find the pseudoinverse of the matrix ```M``` using singular value
decomposition.

> ARGUMENTS
>   - **```M — array_like```**: An ```(m x n)``` matrix whose pseudoinverse we
>     want to find.
>   - **```return_condition_number — Boolean```**: If set to ```True```,
>     function will return the condition number of the pseudoinverse as well.
>     Defaults to ```False```.
>
> OUTPUT
>   - **```pinv```**: An ```(n x m)``` matrix, the pseudoinverse of ```M```.
>   - **```condition_number```**: The condition number of the pseudoinverse.
>     This is only returned if ```return_condition_number``` is set to
>     ```True```.


<a name='vessel'></a>
### ```utils.Vessel```

A convenient container object for holding data, including lists, arrays,
dictionaries, etc. The object can be easily saved and loaded.


> **METHODS**
   
> **```Vessel(filename=None)```**

> Creates a new ```Vessel``` object.   

> > ARGUMENTS
> > - **```filename — str```**: An optional filename for where the data should
> >   be saved. If a filename is specified, the Vessel object will look for a
> >   file at that location; if a valid file is found, it will load the
> >   contents of that file into the object. Otherwise, it will set the
> >   object's current_filename to filename.
>
> > OUTPUT    
> > - None 

> **```Vessel.save(filename=None)```**

> Saves the content of the object to the specified filename using cPickle's
> highest protocol. If no filename is specified, the object will use its 
> current filename. If none exists, an error is thrown.

> > ARGUMENTS
> > - **```filename — str```**: An optional filename for where the data should
> >   be saved. 
>
> > OUTPUT    
> > - None 
>
> **```Vessel.load(filename=None)```**

> Loads the content of the object in the specified using cPickle's
> highest protocol. If no filename is specified, the object will use its 
> current filename. If none exists, an error is thrown.

> > ARGUMENTS
> > - **```filename — str```**: An optional filename for where the data should
> >   be saved. 
>
> > OUTPUT    
> > - None 
>
> **PROPERTIES** <a name='fit-properties'></a>
> 
> When the ```Vessel``` class is instantiated, it makes several properties
> available to the user.
>
>   - **```keys - list```**: The currently available properties on the object.
>   - **```current_filename - str```**: The current filename of the object. If
>     the save method is called without input, this filename will be used.
>
> **USAGE** <a name='vessel-usage'></a>
>
> The Vessel class is a convenient way to store, save, and reload data. Under
> the hood it uses cPickle to efficiently store a wide range of data types. To
> get a sense of what it can do, check out the following examples.
>
> ```python
> import numpy as np
> from mathtools.utils import Vessel
> 
> # Generate some noisy data to fit.
> experiment_name = 'sine wave'
> x = np.linspace(0, 3*np.pi, 200)
> y = np.cos(x/4) + np.random.randn(len(x))
> 
> # Create a vessel to store this data.
> v = Vessel('my_sine_data.dat') # note that extension is arbitrary.
> v.name = experiment_name
> v.x = x
> v.y = y
> 
> # Save the data.
> v.save()
>
> # Many years pass...
>
> g = Vessel('my_sine_data.dat')
> g.name # => 'sine wave'
> g.x    # ==> array([ 0.        ,  0.04736069,  0.09472139, ...
> # etc...
