"""A setuptools based setup module.
See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

# load version from source code
from vireadb import VERSION

here = path.abspath(path.dirname(__file__))
setup(
    name='vireadb',  # Required
    version=VERSION,  # Required
    description='vireadb: Viral Read Database',  # Required
    long_description='ViReaDB is a user-friendly database for storing reference-compressed viral sequence data and computing consensus genome sequences.',  # Optional
    long_description_content_type='text/plain',  # Optional (see note above)
    url='https://github.com/niemasd/vireadb',  # Optional
    author='Niema Moshiri',  # Optional
    author_email='niemamoshiri@gmail.com',  # Optional
    classifiers=[  # Optional
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',

        # Pick your license as you wish
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    keywords='virus sequence database',  # Optional
    packages=find_packages(exclude=['contrib', 'docs', 'tests']),  # Required
    extras_require={  # Optional
        'dev': ['check-manifest'],
    },
    project_urls={  # Optional
        'Bug Reports': 'https://github.com/niemasd/vireadb/issues',
        'Source': 'https://github.com/niemasd/vireadb',
    },
)
