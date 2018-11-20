from setuptools import setup, find_packages

setup(
    name='bac',

    version='0.0.1.dev1',

    description='Binding affinity calculator',

    long_description='Copy from README file',

    url='http://ccs.chem.ucl.ac.uk',

    author='CCS',

    install_requires=['PyYAML', 'parmed', 'numpy', 'supproperty',
                      'pandas', 'mdtraj', 'freesasa'],

    packages=find_packages(),

    include_package_data=True,
)
