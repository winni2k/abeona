#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from glob import glob
from os.path import splitext, basename

from setuptools import find_packages, setup

with open('README.rst', 'r', encoding='utf-8') as f:
    readme = f.read()

setup(
    name='abeona',
    version='0.19.0',
    description='',
    long_description=readme,
    author='Warren W. Kretzschmar',
    author_email='warrenk@kth.se',
    maintainer='Warren W. Kretzschmar',
    maintainer_email='warrenk@kth.se',
    url='https://github.com/winni2k/abeona',
    license='Apache-2.0',

    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: Implementation :: CPython',
        'Programming Language :: Python :: Implementation :: PyPy',
    ],
    entry_points={
        'console_scripts': ['abeona=abeona.__main__:main'],
    },

    install_requires="""
    cortexpy >= 0.38.0
    snakemake
    pandas
    progressbar2
    """.split('\n'),

    packages=find_packages('src'),
    package_dir={'': 'src'},
    py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],
    include_package_data=True,
    zip_safe=False,
)
