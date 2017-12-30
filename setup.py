#!/usr/bin/env python

from distutils.core import setup

setup(
    name='Maybrain', 
    version='0.5.1',
    author='Martyn Rittman and Timothy Rittman',
    author_email='mrittman@physics.org',
    packages=['maybrain', 'maybrain.algorithms', 'maybrain.utils', 'maybrain.plotting'],
    url='https://github.com/RittmanResearch/maybrain/releases',
    license='LICENSE',
    description='A module written in Python for visualising brain connectivity data. Its purpose is to allow easy '
                'visualisation of brain connectome and related data, and perform various analyses.',
    long_description=open('README.md').read(),
    requires=['networkx', 'numpy']

)
