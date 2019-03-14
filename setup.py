#!/usr/bin/env python3

import sys
import setuptools
import mpf_version

with open("README.md", "r") as fd:
    long_description = fd.read()

setuptools.setup(
    name="pympf",
    version=mpf_version.version,
    author="Florian Schanda",
    author_email="florian@schanda.org.uk",
    description="An arbitrary precision IEEE-754 implementation in Python",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/florianschanda/PyMPF",
    project_urls={
        "Bug Tracker"   : "https://github.com/florianschanda/PyMPF/issues",
        "Documentation" : "https://florianschanda.github.io/PyMPF/",
        "Source Code"   : "https://github.com/florianschanda/PyMPF",
    },
    packages=setuptools.find_packages(),
    python_requires=">=3.5, <4",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Software Development :: Libraries",
        "Topic :: System :: Emulators"
    ],
)
