#!/usr/bin/env python

from distutils.core import setup

setup(
    name='Maybrain', 
    version='0.3.0',
    author='Martyn Rittman and Timothy Rittman',
    author_email='mrittman@physics.org',
    packages=['maybrain', 'maybrain.gui'],
#    scripts=['brainObj.py','extraFns.py', 'highlightObj.py', 'mbplot.py', 'writeFns.py'],
    url='https://github.com/RittmanResearch/maybrain/releases',
    license='LICENSE',
    description='A package for visualizing brain network data',
    long_description=open('README.md').read(),
#    install_requires=[
#        'networkX == v1.8.1',
#        'nibabel >= v1.3.0',
#        'numpy >= v1.8.1',
#        'mayavi == v4.4.2'
#   ],
)

