#To build:
# >>> python setup.py --quiet build_ext --inplace clean --all

import os
import sys
from glob import glob
import setuptools
from setuptools import setup

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.join(DIR, "extern", "pybind11"))

from pybind11.setup_helpers import Pybind11Extension, build_ext

del sys.path[-1]

with open("README.md", "r") as README:
    long_description = README.read()

FC = Pybind11Extension(
    'FC',
    sources=['src/FeldmanCousins/FC.cc'],
    language = 'c++',
    extra_compile_args = ['-std=c++11']
)

ProtonDecayTools = Pybind11Extension(
    'ProtonDecayTools',
    sources=['src/ProtonDecay/BayesianCpp/ProtonDecayTools.cc'],
    language = 'c++',
    extra_compile_args = ['-std=c++11', '-I src/ProtonDecay/BayesianCpp/']
)

setup(
    name="Zstats",
    version="2.0.0",
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
    install_requires=["numpy>=1.19.1", "scipy>=1.4.1","mpmath>=1.1.0"],
    cmdclass={"build_ext": build_ext},
    ext_modules=[FC, ProtonDecayTools]
)
