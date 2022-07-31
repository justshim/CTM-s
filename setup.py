#!/usr/bin/env python

from setuptools import setup

setup(name='CTM-s',
      version='1.0',
      # list folders, not files
      packages=['CTM-s.model',  'CTM-s.model.optimization'],
      scripts=['CTM-s/bin/top_level_script.py'],
      package_data={'CTM-s': ['CTM-s/data']},
      )