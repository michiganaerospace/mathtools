# Math Tools

The ```mathtools```  module provides algorithms for processing and
manipulating data. At the current time, Math Tools contains only the ```fit```
module, which provides algorithms for efficient fitting data using regularized
least squares. For more details on the individual modules, please see below.

### Available modules.

- [Fit](#fit) — The ```fit``` module provides algorithms for reguarlized least
squares using different bases.


## Fit 

The ```fit``` module provides the ```Fit``` object — a machine that allows you
to easily fit data using regularized least squares. The goal is to fuse
flexibility and computational efficiency with an easy-to-use API. 

```python
Fit(x=None, y=None, nb_orders=0, basis_type='legendre', reg_coefs=[0.0, 0.0, 0.0]) 
```

#### Inputs

```x``` — array_like [default: None]
    Vector of abscissa values — an ```nb_samples``` length vector of 'x values'.

```y``` — array_like [default: None]
    Vector of ordinate values — an ```nb_samples``` length vector 'y values'.

```nb_orders``` - int [default: 0]
    The number of basis vectors to use when fitting the data. In the case of a
    cubic spline basis, this corresponds to the number of knots used.

```basis_type``` - string [default: ```'legendre'```]
    The type of basis to use for the fits. May have values of 'legendre', 
    'fourier', or 'cubic-spline'.

#### Outputs

```fit``` — Fit Object
   An instance of a ```Fit``` object.


### Examples

To get a sense of what
```Fit``` can do, consider the following example: We have data ```y``` sampled
at points ```x```, as illustrated in the following figure.

