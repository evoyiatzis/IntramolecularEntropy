#!/usr/bin/env python

from distutils.core import setup

try:
    import numpy
except ImportError:
    raise ImportError('numpy cannot be automatically installed by this ' \
                      'tool. Please, install it manually.')

classifiers = [
    'Development Status :: 5 - Production/Stable',
    'Intended Audience :: Education',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
    'Operating System :: MacOS :: MacOS X',
    'Operating System :: POSIX :: Linux',
    'Programming Language :: Python :: 2',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.6',
    'Topic :: Scientific/Engineering :: Physics'
]

def do_setup():
    setup(name='IntramolecularEntropy',
        version='1.0',
        author='Evangelos Voyiatzis',
        author_email='e.voyiatzis@theo.chemie.tu-darmstadt.de',
        url='https://github.com/evoyiatzis/IntramolecularEntropy',
        description='A tool to compute the intramolecular entropy of a particle' \
                  ' system using the quasi harmonic approximation',
        long_description='  ' ,
        download_url='https://github.com/evoyiatzis/IntramolecularEntropy',
        license='LGPLv3',
        packages=['Source'],
        install_requires=['numpy'],
        classifiers = classifiers,
     )

if __name__ == "__main__":
    do_setup()
