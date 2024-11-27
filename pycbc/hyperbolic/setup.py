"""
setup.py file for hyperbolic waveforms
"""

from setuptools import Extension, setup, Command
from setuptools import find_packages

VERSION = '0.0.dev0'

setup (
    name = 'pycbc-hyperbolic',
    version = VERSION,
    description = 'Hyperbolic Encounters of Compact objects',
    long_description = open('descr.rst').read(),
    author = 'Shrey Aggarwal',
    author_email = 'aggar104@umn.edu',
    url = 'http://www.pycbc.org/',
    download_url = 'https://github.com/aggar104/pycbc/blob/master/pycbc/waveform/hyperbolic',
    keywords = ['pycbc', 'signal processing', 'gravitational waves', 'hyperbolic encounters'],
    install_requires = ['pycbc'],
    py_modules = ['hyperbolic'],
    entry_points = {"pycbc.waveform.td":"revchirp = hyperbolic_td",
                    "pycbc.waveform.fd":"revchirp = hyperbolic_fd"},
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.9',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: NA :: NA',
    ],
)