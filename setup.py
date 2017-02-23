#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup
from setuptools.command.test import test as TestCommand

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'Click>=6.0',
    # TODO: put package requirements here
]

test_requirements = [
    # TODO: put package test requirements here
]


class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = ['-v']
        self.test_suite = True

    def run_tests(self):
        # import here, cause outside the eggs aren't loaded
        import pytest
        pytest.main(self.test_args)


setup(
    name='bitk3',
    version='0.1.0',
    description="Set of functions and scripts for bioinformatics",
    long_description=readme + '\n\n' + history,
    author="Davi Ortega",
    author_email='ortegad@caltech.edu',
    url='https://github.com/daviortega/bitk3',
    scripts=['bin/*bk3'],
    packages=[
        'bitk3',
    ],
    package_data={
        'sampledata': ['*']
    },
    package_dir={'bitk3':
                 'bitk3'},
    entry_points={
        'console_scripts': [
            'bitk3=bitk3.cli:main'
        ]
    },
    include_package_data=True,
    install_requires=requirements,
    license="MIT license",
    zip_safe=False,
    keywords='bitk3',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    # test_suite='tests',
    tests_require=['pytest'],
    cmdclass={'test': PyTest}
)
