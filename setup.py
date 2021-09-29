from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
import glob
import os
import pkg_resources

from beastgenpy import __version__, _program

setup(name='beastgenpy',
      version=__version__,
      packages=find_packages(),
      scripts=[
            "beastgenpy/scripts/core_funcs.py",
            "beastgenpy/scripts/glm_funcs.py",
            "beastgenpy/scripts/taxon_set_funcs.py",
            "beastgenpy/scripts/trait_analysis_funcs.py"
            ],
      install_requires=[
            "mako>=1.1",
        ],
      description='Python command line tool to generate BEAST XMLs',
      url='https://github.com/ViralVerity/beastgenpy',
      author='Verity Hill',
      author_email='verity.hill@ed.ac.uk',
      entry_points="""
      [console_scripts]
      {program} = beastgenpy.command:main
      """.format(program = _program),
      include_package_data=True,
      keywords=[],
      zip_safe=False)
