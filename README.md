=========
Maybrain
=========

Overview
========
Maybrain is a module for visualizing brain connectivity data. It has been written by Timothy Rittman
and Martyn Rittman and we are happy to have others contribute code and feedback to the project.

Installation
============
Maybrain has been tested with Python 2.7 and requires the following modules

* networkX == v1.8.1
* nibabel >= v1.3.0
* numpy >= v1.8.1
* mayavi == v4.4.2

We currently have no idea how it will behave with Python 3.x, but will have to make the
switch at some point.

To install, simply copy the files into a convenient location such as Python[xx]\Lib\site-packages
where [xx] is the python version you are using. We will also soon look at uploadin to Pypi for 
easier access.


Some practical details
--------------------------
Modules are imported into the initialistion file to avoid having to call the 
submodule each time.

To avoid having to import mayavi you can use, e.g.

	from mayavin import brainObj

General layout of the module
=============================

There are two principal objects: one for what we call 'the brain', a network with associated
properties. The other is for plotting. Brains can be highlighted so that a subnetwork can be
worked with instead of the entire object.

Metrics and manipulations are included in the brain module and a few extra functions.

For detailed information, please see the pdf user guide.



Rest format

Some text intro

	some code


*italics* and **bold** and ``monospace``

A section 
=========

A sub-section
-------------

Numbered lists 

1. item here

2. item here


some more text


Lists

* number one

* number two indent
  properly

