
''' srbogrid setup script'''

from setuptools import setup


def readme():
    '''Return the contents of the README.md file.'''
    with open('README.md') as freadme:
        return freadme.read()


setup(
    author="Lukasz Mentel",
    author_email="lmmentel@gmail.com",
    description="Package for calculating space-reduced bond-order grids for diatomics",
    include_package_data=False,
    license=open('LICENSE.rst').read(),
    long_description=readme(),
    name='srbogrid',
    packages=['srbogrid'],
    url='https://github.com/lmmentel/srbogrid',
    version='0.1.0',
    classifiers=[
        'Environment :: Console',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
    ],
)
