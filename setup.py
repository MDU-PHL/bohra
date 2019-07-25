from setuptools import setup, find_packages
import bohra


with open("README.md", "r") as fh:
    long_description = fh.read()

VERSION = bohra.__version__

setup (name = 'bohra',
       version = VERSION,
       classifiers = [
         "Programming Language :: Python :: 3.6",
         "Programming Language :: Python :: 3.7",
         "Operating System :: OS Independent",
         "Development Status :: 4 - Beta ",
         "Intended Audience :: Science/Research",
         "Topic :: Scientific/Engineering :: Bio-Informatics",
         "License :: OSI Approved :: GNU General Public License v3 (GPLv3)"
         ],
       include_package_data = True,
       packages=['bohra'],
       description = 'A bioinformatics pipeline for analysing short read Illumina data microbiological public health.',
       long_description = long_description,
       long_description_content_type="text/markdown",
       author = 'Kristy Horan',
       url = 'https://github.com/kristyhoran/bohra',
       install_requires = ['datetime', 'pytest', 'jinja2','biopython>=1.70','pandas>=0.23.0', 'numpy', 'svgwrite', 'psutil','sh','packaging', 'snakemake>=5.4.0'],
       python_requires='>=3.6',
       entry_points={
        "console_scripts": [
           "bohra=bohra.bohra:main"
        ],},
        setup_requires=['pytest-runner'],
        tests_require = ['pytest'],
        test_suite = 'test'
)
