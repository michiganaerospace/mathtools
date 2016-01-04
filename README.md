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


### Fit 

The ```fit``` module provides the ```Fit``` object — a machine that allows you
to easily fit data using regularized least squares. 


```python
Fit(x=None, y=None, nb_orders=0, basis_type='legendre', reg_coefs=[0.0, 0.0, 0.0]) 
```

#### Arguments

```x — array_like [default: None]```
    
Vector of abscissa values — an ```nb_samples``` length vector of 'x values'.

```y — array_like [default: None]```

Vector of ordinate values — an ```nb_samples``` length vector 'y values'.

```nb_orders — int [default: 0]```

The number of basis vectors to use when fitting the data. In the case of a
cubic spline basis, this corresponds to the number of knots used.

```basis_type — string [default: 'legendre']```

The type of basis to use for the fitting. May have values of ```legendre```,
```fourier```, or```cubic-spline```.

#### Methods

```python
Fit.fit(x=None, y=None)
```

The method fits the data using the current configurations of the ```Fit```
object.

```python
Fit.resample(x=None)
```

Resamples the current fit to the specified domain ```x``` using the existing
coefficients. If any part of ```x``` lies outside the original domain of the
fit, the values for the fit in that region are set to zero.

#### Examples

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
to do this in just one line:

```python
# Fit the noisy data.
f = Fit(t, y, 15)
```

This will create an instance of a ```Fit``` object, generate a Legendre
polynomial basis for the provided abscissa — here, the vector ```t``` — with 15
basis vectors, and fit the data ```y``` using this basis. The results of the
fit are always available in the ```results``` property of the ```Fit```
instance. To see how well our fit worked, let's plot it on top of the noisy
data.

```python
r = f.results
plot(r.x, r.y, linewidth=2)
```

![Noisy Data](https://goo.gl/9ozbXw)

Note that although the blue curve smoothly fits the red data points, we are not
explicitly penalizing the derivatives here. We are instead effectively
regularizing by using a small number of Legendre polynomial basis vectors — in
this case, fifteen basis vectors. For more on this topic, check out the
techical discussion at ...

But what if we wanted to sample this function at a different set of points? As
it stands, our smooth curve is sampled only at the values of ```x```
corresponding to the original data. If we wish to sample our smooth curve on a
different set of data points, we can use the ```Fit``` object's ```resample```
method. For example, suppose we wish to sample on a different number of points,
on a different domain.

```python
t_new = np.linspace(-np.pi, 16*np.pi, 500)
```

Note that we are now sampling at 500 points, extending from -pi to 16pi. Using
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

