#!/usr/bin/env python

from setuptools import setup

version = '1.0.0'

required = open('requirements.txt').read().split('\n')

setup(
    name='sequenceQC_reporter',
    version=version,
    description='sequenceQC_reporter reports QC measurements for single-cell sequencing data and generates a standard report.',
    author='Foad Green',
    author_email='foadgreen@gmail.com',
    url='https://github.com/czbiohub/sequenceQC_reporter',
    packages=['sequenceQC_reporter'],
    install_requires=required,
    long_description='See ' + 'https://github.com/czbiohub/sequenceQC_reporter',
    license='MIT',
    entry_points={"console_scripts": ['sequenceQC_reporter = sequenceQC_reporter.cli:cli']}
)
