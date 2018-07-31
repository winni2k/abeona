#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from glob import glob
from os.path import splitext, basename

from setuptools import find_packages, setup

with open('src/abeona/__init__.py', 'r') as f:
    for line in f:
        if line.startswith('__version__'):
            version = line.strip().split('=')[1].strip(' \'"')
            break
    else:
        version = '0.0.1'

with open('README.rst', 'r', encoding='utf-8') as f:
    readme = f.read()

setup(
    name='abeona',
    version=version,
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
        'console_scripts': ['abeona=abeona.__main__:main_without_argv'],
    },

    install_requires="""
    cortexpy >= 0.31.3
    snakemake
    pandas
    """.split('\n'),

    packages=find_packages('src'),
    package_dir={'': 'src'},
    py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],
    include_package_data=True,
    zip_safe=False,
)
