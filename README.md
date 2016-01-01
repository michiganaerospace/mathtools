# Math Tools

The ```mathtools```  module provides algorithms for processing and
manipulating data. At the current time, Math Tools contains only the ```fit```
module, which provides algorithms for efficient fitting data using regularized
least squares. For more details on the individual modules, please see below.

### Available modules.

- [Fit](#fit) — The ```fit``` module provides algoritms for reguarlized least
squares using different bases.


## Fit 

The ```fit``` module provides the ```Fit``` object, which allows you to easily
fit data using regularized least squares. To understand what it does, consider
the following example. We have data ```y``` sampled at points ```x```.
