from setuptools import setup, find_packages

setup (name = 'bohra',
       version = '1.0.1',
       include_package_data = True,
       packages=find_packages(),
       description = 'A bioinformatics pipeline for analysing short read Illumina data microbiological public health.',
       author = 'Kristy Horan',
       url = 'https://github.com/MDU-PHL/bohra',
       install_requires = ['jinja2','biopython>=1.70','pandas>=0.23.0', 'pathlib', 'numpy', 'snakemake>=5.3.0'],
       python_requires='>=3.6',
       entry_points={
        "console_scripts": [
            "bohra = bohra:main"
        ],},
        setup_requires=['pytest-runner'],
        tests_require = ['pytest'],
        test_suite = 'test'
)
