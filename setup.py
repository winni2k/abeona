import json
import os

from setuptools import find_packages, setup

with open('abeona/__init__.py', 'r') as f:
    for line in f:
        if line.startswith('__version__'):
            version = line.strip().split('=')[1].strip(' \'"')
            break
    else:
        version = '0.0.1'

with open('README.rst', 'r', encoding='utf-8') as f:
    readme = f.read()


def get_requirements_from_pipfile_lock(pipfile_lock=None):
    if pipfile_lock is None:
        pipfile_lock = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'Pipfile.lock')
    lock_data = json.load(open(pipfile_lock))
    return [package_name for package_name in lock_data.get('default', {}).keys()]


REQUIRES = get_requirements_from_pipfile_lock()

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

    keywords=[
        '',
    ],

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

    install_requires=REQUIRES,
    tests_require=['coverage', 'pytest'],

    packages=find_packages(),
)
