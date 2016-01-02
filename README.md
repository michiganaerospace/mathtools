# Math Tools

The ```mathtools```  module provides algorithms for processing and
manipulating data. At the current time, Math Tools contains only the ```fit```
module, which provides algorithms for efficient fitting data using regularized
least squares. For more details on the individual modules, please [see
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


### Available modules

- [Fit](#fit) — The ```fit``` module provides algorithms for reguarlized least
squares using different bases.


## Fit 

The ```fit``` module provides the ```Fit``` object — a machine that allows you
to easily fit data using regularized least squares. 


```python
Fit(x=None, y=None, nb_orders=0, basis_type='legendre', reg_coefs=[0.0, 0.0, 0.0]) 
```

#### Arguments

```x``` — array_like [default: None]
    
Vector of abscissa values — an ```nb_samples``` length vector of 'x values'.

```y``` — array_like [default: None]

Vector of ordinate values — an ```nb_samples``` length vector 'y values'.

```nb_orders``` — int [default: 0]

The number of basis vectors to use when fitting the data. In the case of a
cubic spline basis, this corresponds to the number of knots used.

```basis_type — string [default: 'legendre']```

The type of basis to use for the fits. May have values of ```legendre```,
```fourier```, or```cubic-spline```.

#### Methods

```Fit.fit(x=None, y=None)```

The method fits the data using the current configurations of the ```Fit```
object.

#### Examples

To get a sense of what ```Fit``` can do, consider the following example: We
have data ```y``` sampled at points ```x```, as illustrated in the following
figure.

For example, suppose we have some noisy data, ```y``` sampled at points
```x```.  We wish to fit this data using smoothed Legendre polynomials. Well,
it's pretty easy.

```python
f = Fit(x, y, nb_orders=15, reg_coefs=[0,1e-1, 1e-2])
f.config(x, nb_orders, reg_coefs)
f.fit(x,y)
```
