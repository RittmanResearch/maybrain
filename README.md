# Maybrain


[![Build Status](https://travis-ci.org/RittmanResearch/maybrain.svg?branch=master)](https://travis-ci.org/RittmanResearch/maybrain)[![Python v3.5](https://img.shields.io/badge/python-v3.5-blue.svg)]() [![Release](https://img.shields.io/github/release/rittman/maybrain/all.svg)](https://github.com/rittman/maybrain/releases) [![GitHub Issues](https://img.shields.io/github/issues/rittman/maybrain.svg)](https://github.com/rittman/maybrain/issues) [![License CC-BY-4.0](https://img.shields.io/badge/license-CC--BY--4.0-lightgrey.svg)](https://github.com/rittman/maybrain/blob/master/LICENSE)


Maybrain is a module written in Python for visualising brain connectivity data. Its purpose is to allow easy visualisation of brain connectome and related data, and perform various analyses. 

## Dependencies

To run Maybrain you will need a Python 2.7 installtion and several other packages on which parts of the code depend. The following are required for analysis:

* [Numpy](http://www.numpy.org/) 1.8.1
* [NetworkX](http://networkx.github.io/) 1.8.1

The following is required for plotting functions:
* [Mayavi](http://docs.enthought.com/mayavi/mayavi/) 4.4.2

The following provides some extra functionality for input of certain data types:
* [NiBabel](http://nipy.org/nibabel/) 1.3.0

Software versions higher than those given may also work. If you are not familiar with Python, or you don't want to manually install and deal with each package separately, it is recommended that you install a pre-packaged version that will include most of the above, for example [Anaconda](https://www.anaconda.com) or [Enthought Canopy](https://www.enthought.com/downloads/). Installation for each package can be found on the individual websites.

We currently have no idea how it will behave with Python 3.x, but will have to make the switch at some point.

## Installation

To install, simply copy the files that are inside the `maybrain` folder into the **site-packages** folder of your Python installation.

If you are not sure where this folder is located, you can open a Python terminal and run the following code to find where Python searches for installed packages:

```python
import sys
print('\n'.join(sys.path))
```

## Documentation
For a detailed documentation, with usage examples and explanations, please go to our [Wiki pages](https://github.com/rittman/maybrain/wiki).


## Contributing
The authors are happy for developers to extend, customise, improve or simply to [create an issue](https://github.com/rittman/maybrain/issues). We will create an innovative and meaningful hall of fame for anyone contributing a bug-fix and promise to buy you a beer (or acceptable non-alcoholic alternative) when we meet.

### Submitting a Pull Request
1. Fork it.
2. Create a branch (`git checkout -b my_maybrain`)
3. Commit your changes (`git commit -am "Added message type"`)
4. Push to the branch (`git push origin my_maybrain`)
5. Open a [Pull Request](https://github.com/rittman/maybrain/pulls)
6. Enjoy a refreshing beverage and wait

## Credits

Maybrain has been written by Timothy Rittman and Martyn Rittman.

## License

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
