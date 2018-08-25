#!/usr/bin/env python

from setuptools import setup

setup(
    name='maybrain',
    version='0.7.3',
    author='Martyn Rittman, Timothy Rittman and Tiago Azevedo',
    author_email='mrittman@physics.org',
    packages=['maybrain', 'maybrain.algorithms', 'maybrain.utils', 'maybrain.plotting', 'maybrain.resources'],
    url='https://github.com/RittmanResearch/maybrain',
    license='LICENSE',
    description='A module written in Python for visualising brain connectivity data. Its purpose is to allow easy '
                'visualisation of brain connectome and related data, and perform various analyses.',
    long_description=open('README.md').read(),
    requires=['networkx', 'numpy', 'matplotlib', 'nilearn'],
    include_package_data=True
)
