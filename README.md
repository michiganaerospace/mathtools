# Math Tools

The ```mathtools```  package provides algorithms for processing and manipulating
data. At the current time, it contains only the ```fit``` module, which
provides algorithms for efficiently fitting data using regularized least
squares.  For more details on the individual modules, please [see
below](#available-modules).

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
```/tests``` directory, and can be run using the nose testing framework. For a
good introduction to nose, check out this
[resource](http://goo.gl/hfzYPz). 

As a general rule, before any code is checked in to the repository, one should
verify that all unit tests are passing.

## Tools 

- [fit.Fit](#fitfit) — The ```Fit``` class provides algorithms for regularised
  least squares using different bases.
    - [methods](#methods)
    - [properties](#properties)
    - [usage](#usage)
- [legendre.legendre_basis](#legendre) — Function for generating a Legendre
polynomial basis.
- [legendre.d_legendre_basis](#d_legendre) — Generates the first derivative of 
a Legendre polynomial basis.
- [legendre.d2_legendre_basis](#d2_legendre) — Generates the second derivative
of a Legendre polynomial basis.
- [create_legendre_basis](#create_legendre) — Generate a Legendre basis
  object.

### ```fit.Fit```

The ```Fit``` class is a machine that allows you to easily fit data using
regularized least squares. It provides a convenient interface to the core basis
generation and fitting routines. The details of these routines may be found at
[legendre.create_legendre_basis](#legendrecreatelegendrebasis).

#### Methods

##### ```Fit.__init__(x=None, nb_bases=0, basis_type='legendre', reg_coefs=[0.0, 0.0, 0.0])```

To create a fit object requires, at a minimum, the domain and the number of
bases to be specified.

> ARGUMENTS    
>   - **```x — array_like```**: Vector of abscissa values — an ```nb_samples```
>     length vector of 'x values'.
>   - **```nb_bases — int```**: The number of basis vectors to use when fitting
>     the data.  In the case of a cubic spline basis, this corresponds to the
>     number of knots used.
>   - **```basis_type — str```**: The type of basis to use for the fitting. May
>     have values of ```legendre```, ```fourier```, or```cubic-spline```.
>   - **```reg_coefs — array_like```**: A list or array of three regularization
>     coefficients for penalizing the magnitude of the fit and its first and
>     second derivatives, respectively. The default value is ```reg_coefs=[0.0,
>     0.0, 0.0]```.
>
> OUTPUTS    
>   - None

##### ```Fit.fit(y)```

The method fits the data, ```y```, using the current configurations of the ```Fit```
object.

> ARGUMENTS   
>   - **```y — array_like```**: Vector of data samples that we want to fit. The
>     in data in the vector ```y``` must correspond to the domain ```x``` for
>     which the basis was constructed.
>
> OUTPUTS  
>   - **```results — object```**: An object containing fit results. It has the
>     following properties:
>       - ```x — array_like```: The domain associated with the fit.
>       - ```y — array_like```: The fit sampled on the domain. 
>       - ```dy — array_like```: The derivative of the fit.
>       - ```d2y — array_like```: The second derivative of the fit.
>       - ```coefs — array_like```: The coefficients of the fit.

##### ```Fit.resample(x)```

Resamples the current fit to the specified domain ```x``` using the existing
coefficients. If any part of ```x``` lies outside the original domain of the
fit, the values for the fit in that region are set to zero. The ```fit()``` 
must have been run before the fit can be resampled.

> ARGUMENTS    
>   - **```x — array_like```**: A vector of abscissa values — an
>     ```nb_samples``` length vector of x-values.
>
> OUTPUTS  
>   - **```results — object```**: An object containing fit results. It has the
>     following properties:
>       - ```x — array_like```: The domain associated with the fit.
>       - ```y — array_like```: The fit sampled on the domain.
>       - ```dy — array_like```: The derivative of the fit.
>       - ```d2y — array_like```: The second derivative of the fit.
>       - ```coefs — array_like```: The coefficients of the fit.

##### ```Fit.config(x=None, nb_bases=0, basis_type=None, reg_coefs=None)```

Recomputes the basis given a change in the underlying parameters. When the 
basis is update, all existing fit results are discarded. 

> ARGUMENTS    
>   - **```x — array_like```**: Vector of abscissa values — an ```nb_samples```
>   - **```nb_bases — int```**: The number of basis vectors to use when fitting
>     the data.  In the case of a cubic spline basis, this corresponds to the
>     number of knots used.
>   - **```basis_type — str```**: The type of basis to use for the fitting. May
>     have values of ```legendre```, ```fourier```, or```cubic-spline```.
>   - **```reg_coefs — array_like```**: A list or array of three regularization
>     coefficients for penalizing the magnitude of the fit and its first and
>     second derivatives, respectively. The default value is ```reg_coefs=[0.0,
>     0.0, 0.0]```.
>
> OUTPUTS   
>   - None.

#### Properties

When the ```Fit``` class is instantiated, it makes several properties available
to the user.

>   - **```basis — object```**: The current basis object. For details on what
>     basis objects are and how to use them, check out the discussion of [basis
>     objects](#basis_object).
>   - **```basis_type — str```**: The type of basis to use for the fitting. May
>     have values of ```legendre```, ```fourier```, or```cubic-spline```.
>   - **```coefs — array_like```**: The current fit coefficients. If there is
>     no active fit, the coefficients are set to  ``None``.
>   - **```nb_bases — int```**: The number of basis vectors currently used in
>     the basis.
>   - **```reg_coefs — array_like```**: A list or array of three regularization
>   coefficients for penalizing the magnitude of the fit and its first and
>   second derivatives, respectively. 
>   - **```x — array_like```**: The domain on which the basis was built.



#### Usage

To get a sense of what the ```Fit``` class can do, let's try to fit some noisy
data. We'll generate a sine wave and add a little noise to it, as illustrated
in the following code snippet.

```python
# Create some noisy data.
t = np.linspace(0,15*np.pi, 300)
y = np.sin(2*np.pi/5*t) + 0.2 * np.random.randn(len(t))
```

The noisy sine data is shown in the following figure.

![Noisy Data](https://goo.gl/elq37W)

We'd like to fit this data with a smooth curve. We can use the ```Fit``` class
to do this quite easily. We first create an instance of the ```Fit``` class,
intializing it with the domain on which we'd like to fit the data, and the
number of basis vectors we want to use, 

```python
# Create a fit object with 15 basis vectors. 
f = Fit(t, 15)
```

This will create an instance of a ```Fit``` object, generate a Legendre
polynomial basis for the provided abscissa (the Legendre polynomial basis is
the default; others may be specified) — here, the vector ```t``` — with 15
basis vectors. To fit the data, we use the object's ```fit``` method,

```python
# Fit the noisy data.
r = f.fit(y)
```

The fit method returns a structure, here ```r```, which contains the fit 
coefficients, the fit itself, and its first two derivatives. For convenience,
the original domain is also included in the structure. To see how well we did,
let's plot the fit on top of the original, noisy data.

```python
plot(r.x, r.y, linewidth=2)
```

![A Good Fit](https://goo.gl/9ozbXw)

Note that although the blue curve smoothly fits the red data points, we are not
explicitly penalizing the derivatives here (the regularization coefficients
are, by default, zero). We are instead *effectively* regularizing by using a
small number of Legendre polynomial basis vectors — in this case, fifteen basis
vectors. For more on this topic, check out the techical discussion at ...

But what if we wanted to sample this function at a different set of points? As
it stands, our smooth curve is sampled only at the values of ```x```
corresponding to the original data. If we wish to sample our smooth curve on a
different set of data points, we can use the ```Fit``` object's ```resample```
method. For example, suppose we wish to sample on a larger number of points,
and on a different domain.

```python
t_new = np.linspace(-np.pi, 6*np.pi, 500)
```

Note that we are now sampling at 500 points, extending from -pi to 6pi. Using
the previously computed ```Fit``` object, we compute a new results object,

```python 
rs = f.resample(t_new)
```

The new fit now samples the same curve (i.e., using the same fit coefficients),
but at different points. In this example, we have deliberately resampled the
curve on a domain that extends beyond the support of the original data. Of
course, we cannot expect to fit data in this region (i.e. the ```Fit``` object
does not extrapolate). The convention here is to set the fit to zero in those
regions that do not intersect the original domain. We illustrate this by
plotting the resampled data:

```python
plot(rs.x, rs.y, linewidth=2)
```

![Resampled Data](https://goo.gl/uxq5ju)

Where the resampled domain intersects the support of the original data, we
reproduce the fit. However, once we venture beyond the support of that data,
the fit returns zero.

As mentioned above, the ```results``` object returned from the ```fit``` method
also contains the derivatives of the fit, which may be useful. To take a look
at the first derivative, for example, we can plot ```r.dy```:

```python
plot(r.x, r.dy, linewidth=2)
```

![Derivative of Fit](https://goo.gl/t5Wjms)

We may want to alter our fit in some way — perhaps by choosing a different 
basis, or selecting a different number of basis vectors. The ```Fit``` class's
 ```config``` method allows us to easily change parameters and recompute the 
basis. For example, suppose we wished to use only ten basis vectors, rather
than fifteen, as above. We simply call ```config()``` with the parameters we 
wish to change:

```python
f.nb_bases # ==> 15
f.config(nb_bases=10)
f.nb_bases # ==> 10
```

Other parameters can be changed similarly. When ```config``` is called, the
basis is recomputed and all fit coefficients are discarded.

<a name="legendre"></a>
### ```legendre.legendre_basis(x, nb_bases)```

Computes the Legendre polynomial basis on a specified domain.

> ARGUMENTS    
>   - **```x — array_like```**: The domain over which we are defining the
>     basis. An ```nb_samples``` length vector.
>   - **```nb_bases — int```**: The number of basis vectors to generate.
>
> OUTPUTS   
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
> OUTPUTS   
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
> OUTPUTS   
> - **```d2B — array_like```**: An ```nb_samples x nb_bases``` array. The
>   columns of the array are the second derivative of the Legendre polynomial
>   basis vectors.


<a name='create_legendre'></a>
### ```legendre.create_legendre_basis(x, nb_bases, reg_coefs=[0,0,0], x_ref=None)```

This function creates a Legendre polynomial basis object. The basis object 
contains everything one needs to perform a regularized least squares fit of
data on the specified domain.

> ARGUMENTS    
>   - **```x — array_like```**: The domain over which we are defining the
>     basis. An ```nb_samples``` length vector.
>   - **```nb_bases — int```**: The number of basis vectors to generate.
>   - **```reg_coefs — array_like```**: A list or array of three regularization
>     coefficients for penalizing the magnitude of the fit and its first and
>     second derivatives, respectively. 
>   - **```x_ref — array_like```**: An optional reference domain. This is
>     useful for resampling data.  It ensures that data is mapped to the
>     interval [-1, 1] in a consistent same way, and allows us to avoid
>     attempting to fit data outside of the original domain.
>
> OUTPUTS <a name='basis_object'></a>   
>   - **```basis — object```**: A basis object. Basis objects consolidate all
>     the information required to fit data with the basis. The basis object
>     contains the following properties and methods:
>       - **```B — array_like```**: An ```nb_samples x nb_bases``` array of
>         Legendre polynomial basis (column) vectors.
>       - **```B_ — array_like```** The so-called 'brick'. The brick is a
>         concatenation of the ```B```, ```I```, ```dB```, ```d2B``` matrices.
>         The brick matrix allows us to force the solution and its derivative
>         to be close to zero, as dictated by the regularization coefficients.
>         In fact, only those components with nonzero regularization
>         coefficients
>         are included in the brick, in order to minimize computational
>         overhead.
>       - **```dB — array_like```**: Derivative of basis vectors in ```B```.
>         ```nb_samples x nb_bases``` in size.
>       - **```d2B — array_like```** : Second derivative of basis vectors in
>         ```B```.   
>       - **test**
