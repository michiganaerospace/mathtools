# Math Tools

The ```mathtools```  module provides algorithms for processing and
manipulating data. At the current time, Math Tools contains only the ```fit```
module, which provides algorithms for efficient fitting data using regularized
least squares. For more details on the individual modules, please see below.

### Available modules.

- [Fit](#fit) — The ```fit``` module provides algorithms for reguarlized least
squares using different bases.


## Fit 

The ```fit``` module provides the ```Fit``` object, — a machine, which allows
you to easily fit data using regularized least squares. The goal is to fuse
flexibility and computational efficiency with an easy-to-use API. 

```python
Fit(x=None, y=None, nb_orders=0, basis_type='legendre',  
    reg_coefs=[0.0, 0.0, 0.0], existing_basis=None, filename=None):
```

This returns an instance of a ```Fit``` object. The fit object allows you to 
easily fit data using regularized least squares.

### Examples

To get a sense of what
```Fit``` can do, consider the following example: We have data ```y``` sampled
at points ```x```, as illustrated in the following figure.

