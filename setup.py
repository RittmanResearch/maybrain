# -*- coding: utf-8 -*-

# to create MANIFEST file, python setup.py sdist

from distutils.core import setup

setup(
    name='Maybrain', 
    version='0.2.2',
    author='Martyn Rittman and Timothy Rittman',
    author_email='mrittman@physics.org',
    packages=['maybrain', 'maybrain.test'],
    scripts=['brainObjs.py','plot.py', 'extraFns.py','recipes.py','writeFns.py'],
    url='http://code.google.com/p/maybrain/',
    license='LICENSE.txt',
    description='A package for visualizing brain network data',
    long_description=open('README.txt').read(),
    install_requires=[
        'networkX == v1.8.1',
        'nibabel >= v1.3.0',
        'numpy >= v1.8.1',
        'mayavi == v4.4.2'
    ],
)

