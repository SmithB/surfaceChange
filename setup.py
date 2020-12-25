import os
from setuptools import setup, find_packages

# get long_description from README.md
with open("README.md", "r") as fh:
    long_description = fh.read()


# list of all scripts to be included with package
scripts = [os.path.join('scripts',f) for f in os.listdir('scripts')] +\
    [os.path.join('surfaceChange', f) for f in ['ATL11_to_ATL15.py']]

setup(
    name='surfaceChange',
    version='0.0.0.1',
    description='Driver for elevation-change mapping in Python.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/bsmith/surfaceChange',
    author='Ben Smith',
    author_email='besmith@uw.edu',
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
    ],
    keywords='altimetry, least-squares, yak_shaving',
    packages=find_packages(),
    scripts=scripts,
    include_package_data=True,
    package_data={'surfaceChange':['resources/*']}
)
