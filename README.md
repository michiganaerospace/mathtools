# Math Tools

The ```mathtools```  module provides algorithms for processing and manipulating
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

f = Fit(x,y,nb_orders=15)
```

## Tests

All code in the ```mathtools``` is fully tested. The tests are available in the
```/tests``` directory, and can be run using the nose testing framework. For a
good introduction to nose, check out this
[resource](http://pythontesting.net/framework/nose/nose-introduction/). 

As a general rule, before any code is checked in to the repository, one should
verify that all unit tests are passing.

## Available modules

- [Fit](#fit) — The ```fit``` module provides algorithms for reguarlized least
squares using different bases.
    - [Methods](#methods)
    - [Properties](#properties)


### Fit 

The ```fit``` module provides the ```Fit``` object — a machine that allows you
to easily fit data using regularized least squares. 

#### Methods

##### ```Fit.__init__(x=None, nb_orders=0, basis_type='legendre', reg_coefs=[0.0, 0.0, 0.0])```
##### Arguments

*```x — array_like```*

Vector of abscissa values — an ```nb_samples``` length vector of 'x values'.

###### ```nb_bases — int```

The number of basis vectors to use when fitting the data. In the case of a
cubic spline basis, this corresponds to the number of knots used.

###### ```basis_type — str```

The type of basis to use for the fitting. May have values of ```legendre```,
```fourier```, or```cubic-spline```.

###### ```reg_coefs — array_like```

A list or array of three regularization coefficients for penalizing the
magnitude of the fit and its first and second derivatives, respectively. The
default value is ```reg_coefs=[0.0, 0.0, 0.0]```

##### ```Fit.fit(y)```

The method fits the data using the current configurations of the ```Fit```
object.

###### ```y — array_like```

Vector of data samples that should be fit. 

##### ```Fit.resample(x)```

Resamples the current fit to the specified domain ```x``` using the existing
coefficients. If any part of ```x``` lies outside the original domain of the
fit, the values for the fit in that region are set to zero.

##### ```Fit.config(x=None, nb_bases=0, basis_type=None, reg_coefs=None)```

Reconfigure the basis. Change the number of basis vectors, the type, or how
much regularization is used.

###### ```x — array_like```

Vector of abscissa values — an ```nb_samples``` length vector of 'x values'.

###### ```nb_bases — int```

The number of basis vectors to use when fitting the data. In the case of a
cubic spline basis, this corresponds to the number of knots used.

###### ```basis_type — str```

The type of basis to use for the fitting. May have values of ```legendre```,
```fourier```, or```cubic-spline```.

###### ```reg_coefs — array_like```

A list or array of three regularization coefficients for penalizing the
magnitude of the fit and its first and second derivatives, respectively. The
default value is ```reg_coefs=[0.0, 0.0, 0.0]```

### Examples

To get a sense of what ```Fit``` can do, let's try to fit some noisy data. We 
generate a sine wave, and add a little noise to it, as illustrated in the 
following code.

```python
# Create some noisy data.
t = np.linspace(0,15*np.pi, 300)
y = np.sin(2*np.pi/5*t) + 0.2 * np.random.randn(len(t))
```

The noisy sine data is shown in the following figure.

![Noisy Data](https://goo.gl/elq37W)

We'd like to fit this data with a smooth curve. We can use the ```Fit``` object
to do this quite easily. We first create an instance of the fit object,
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

![Noisy Data](https://goo.gl/9ozbXw)

Note that although the blue curve smoothly fits the red data points, we are not
explicitly penalizing the derivatives here (the regularization coefficients
are, by default, zero). We are instead effectively regularizing by using a
small number of Legendre polynomial basis vectors — in this case, fifteen basis
vectors. For more on this topic, check out the techical discussion at ...

But what if we wanted to sample this function at a different set of points? As
it stands, our smooth curve is sampled only at the values of ```x```
corresponding to the original data. If we wish to sample our smooth curve on a
different set of data points, we can use the ```Fit``` object's ```resample```
method. For example, suppose we wish to sample on a different number of points,
on a different domain.

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
does not extrapolate). The convention here is to set the data to zero in the
regions that do not intersect the original domain. We illustrate this by
plotting the resampled data:

```python
plot(rs.x, rs.y, linewidth=2)
```

![Noisy Data](https://goo.gl/uxq5ju)

Where the resampled domain intersects the support of the original data, we
reproduce the fit. However, once we leave the support of that data, the fit 
simply returns the data.

As mentioned above, the object returned from the fit method also contains the
derivatives of the fit, which may be useful. To take a look at the  first
derivative, for example, we can plot ```r.dy```:


```python
plot(r.x, r.dy, linewidth=2)
```

![Derivative of Fit](https://goo.gl/b75E2u)



