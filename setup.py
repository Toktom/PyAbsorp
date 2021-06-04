# -*- coding: utf-8 -*-
"""
PyAbsorp setup file
=================
@Author: Michael Markus Ackermann
"""

from setuptools import setup

settings = {
    'name': 'PyAbsorp',
    'version': '0.2.0',
    'description': 'Sound absorption coefficient models implemented in Python.',
    'url': 'https://github.com/Toktom/PyAbsorp',
    'author': 'Michael Markus Ackermann',
    'author_email': 'dev.toktom@outlook.com',
    'license': 'MIT',
    'install_requires': [
        'numpy>=1.20.3', 
        'scipy>=1.6.3',
        'pydantic>=1.8.2'],
    'packages': ['pyabsorp', 'pyabsorp.models', 'pyabsorp.air']
    }

setup(**settings)
