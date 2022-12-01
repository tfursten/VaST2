from setuptools import setup, find_packages
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="vast2",
    version="2.0.0",
    description='Variant Site Strain Typer',
    install_requires=[
        'Click',
    ],
    packages=['vast'],
    entry_points='''
        [console_scripts]
        vast=vast.main:cli
    ''',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/tfursten/VaST',
    author='Tara Furstenau',
    author_email='tara.furstenau@nau.edu',
    classifiers=[ 
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ], 
    python_requires='>=3.6'
)
