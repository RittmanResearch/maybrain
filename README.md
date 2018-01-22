# Maybrain


[![Python v3.5, v3.6](https://img.shields.io/badge/python-v3.5,v3.6-blue.svg)]() [![Build Status](https://travis-ci.org/RittmanResearch/maybrain.svg?branch=master)](https://travis-ci.org/RittmanResearch/maybrain) [![codecov](https://codecov.io/gh/RittmanResearch/maybrain/branch/master/graph/badge.svg)](https://codecov.io/gh/RittmanResearch/maybrain) [![Release](https://img.shields.io/github/release/RittmanResearch/maybrain/all.svg)](https://github.com/RittmanResearch/maybrain/releases) [![License Apache-2.0](https://img.shields.io/github/license/RittmanResearch/maybrain.svg)](https://github.com/rittman/maybrain/blob/master/LICENSE)


Maybrain is a Python package for analysing and visualising brain connectome and related data. 

## Dependencies

To run Maybrain you will need a Python 3.5 installation and several other packages on which parts of the code depend. The following are required for analysis:

* [Numpy](http://www.numpy.org/) 1.13.1
* [NetworkX](http://networkx.github.io/) 2.0

The following is required for plotting functions:
* [Matplotlib](http://matplotlib.org/) 2.1

The following provides some extra functionality for input of certain data types:
* [NiBabel](http://nipy.org/nibabel/) 2.1.0
* [Mayavi](http://docs.enthought.com/mayavi/mayavi/) 4.5.0
* [bctpy](https://github.com/aestrivex/bctpy) 0.5.0

Software versions higher than those given may also work. If you are not familiar with Python, or you don't want to manually install and deal with each package separately, it is recommended that you install a pre-packaged version that will include most of the above, for example [Anaconda](https://www.anaconda.com). Installation for each package can be found on the individual websites, but if you want to use [Anaconda](https://www.anaconda.com) you will find instructions in our [documentation page](https://github.com/RittmanResearch/maybrain/wiki).


## Installation

To install, simply run `setup.py` throughout `pip`:

```bash
$ pip install .
```

If you want more detailed instructions about installing maybrain using Anaconda, please go to our [Wiki pages](https://github.com/RittmanResearch/maybrain/wiki).

## Documentation
For a detailed documentation, with usage examples and explanations, please go to our [Wiki pages](https://github.com/RittmanResearch/maybrain/wiki).


## Contributing
The authors are happy for developers to extend, customise, improve or simply to [create an issue](https://github.com/RittmanResearch/maybrain/issues). We will create an innovative and meaningful hall of fame for anyone contributing a bug-fix and promise to buy you a beer (or acceptable non-alcoholic alternative) when we meet.

### Submitting a Pull Request
1. Fork it.
2. Create a branch (`git checkout -b my_maybrain`)
3. Commit your changes (`git commit -am "Added message type"`)
4. Push to the branch (`git push origin my_maybrain`)
5. Open a [Pull Request](https://github.com/RittmanResearch/maybrain/pulls)
6. Enjoy a refreshing beverage and wait

## Credits

Maybrain was originally written by Timothy Rittman and Martyn Rittman, but with substantial improvements and updates from Tiago Azevedo.


## Licence

Maybrain is licensed under the Apache License, Version 2.0. See [LICENSE](https://github.com/RittmanResearch/maybrain/blob/master/LICENSE) for the full license text.

