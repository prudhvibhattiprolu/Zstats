#To build:
# >>> python setup.py --quiet build_ext --inplace clean --all

import os
import sys
from glob import glob
import setuptools
from setuptools import setup

with open("README.md", "r") as README:
    long_description = README.read()

setup(
    name="Zstats",
    version="2.0",
    author="Prudhvi Bhattiprolu",
    author_email="prudhvibhattiprolu@gmail.com",
    description="Statistical measures for discovery and exclusion of new physics signals at (multi-channel) counting experiments",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT",
    url="https://github.com/prudhvibhattiprolu",
    packages=setuptools.find_packages(where="src"),
    package_dir={"": "src"},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
    install_requires=["numpy>=1.19.1", "scipy>=1.4.1","mpmath>=1.1.0"]
)
