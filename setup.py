from setuptools import setup
from setuptools import find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='polcalibration',
      version='0.0.1',
      url='https://github.com/miguelcarcamov/polcal_scripts',
      description='Object oriented scripts for polarization calibration',
      author='Miguel Carcamo',
      author_email='miguel.carcamo@manchester.ac.uk',
      long_description=long_description,
      long_description_content_type="text/markdown",
      license='GNU',
      packages=find_packages(),
      classifiers=[
          "Programming Language :: Python :: 3",
          "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
          "Operating System :: OS Independent"],
      )
