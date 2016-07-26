#!/usr/bin/env python

from distutils.core import setup

setup(name='OVarCall',
      version='0.1.0',
      description='Python tools to identify somatic mutations, using overlapping read information.',
      author='Takuya Moriyama',
      author_email='moriyama@hgc.jp',
      url='',
      package_dir={'': 'lib'},
      packages=['OVarCall', 'OVarCall.filter', 'OVarCall.methods', 'OVarCall.utilOVar'],
      scripts=['OVarCall'],
      license='GPL-3'
      )
