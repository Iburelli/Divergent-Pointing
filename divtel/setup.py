#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
# import sys
from setuptools import setup, find_packages
import os
import re


def get_property(prop, project):
    result = re.search(r'{}\s*=\s*[\'"]([^\'"]*)[\'"]'.format(prop), open(project + '/__init__.py').read())
    return result.group(1)


setup(name='divtel',
      version=0.1,
      description="Divergent pointing library",  # these should be minimum list of what is needed to run
      packages=find_packages(),
      install_requires=['astropy',
                        ],
      package_data={},
      tests_require=['pytest', 'pytest-ordering'],
      author='T. Vuillaume, A. Donini, T. Gasparetto',
      author_email='thomas.vuillaume[at]lapp.in2p3.fr',
      license='MIT',
      url='',
      long_description='',
      )
