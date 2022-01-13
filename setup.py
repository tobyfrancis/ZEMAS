from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='zemas',
    version='0.0.1',    
    description='ZoNexus Electron Microscopy Acquisition Suite',
    url='https://github.com/tobyfrancisv/ZEMAS.git',
    author='Toby Francis',
    author_email='tobyfrancisv@gmail.com',
    license='GPU',
    packages=find_packages(exclude=['data']),
    install_requires=['PyQt6',
                      'monty',
                      'pymatgen'                     
                      ],

    classifiers = [
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: GNU General Public License v3.0 (GPL3)',       
        'Programming Language :: Python :: 3',
    ],
    long_description = long_description
)