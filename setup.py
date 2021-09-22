from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
import glob
import os
import pkg_resources


setup(name='beastgenpy',
      version="0.1",
      packages=find_packages(),
      scripts=[
            "beastgenpy/scripts/beastgen.py",
            "beastgenpy/scripts/core_funcs.py",
            "beastgenpy/scripts/glm_funcs.py",
            "beastgenpy/scripts/taxon_set_funcs.py",
            "beastgenpy/scripts/trait_analysis_functions.py"
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
      """.format(program = "beastgenpy"),
      include_package_data=True,
      keywords=[],
      zip_safe=False)
