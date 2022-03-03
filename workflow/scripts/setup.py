import setuptools

setuptools.setup(
    name='genetools',
    version="0.0.1",
    description='Download representative/reference genomes from NCBI RefSeq database then filter and compile specific genes from them',
    author='Charlie Bushman',
    author_email='ctbushman@gmail.com',
    url='https://github.com/PennChopMicrobiomeProgram',
    packages=['genetools'],
    entry_points={
        'console_scripts': [
            'genetools=genetools.command:main',
        ],
    },
    install_requires=[
        'tqdm',
        'wget',
    ],
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',
        'Operating System :: POSIX :: Linux',
        'Operating System :: MacOS :: MacOS X',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    license='GPLv2+',
)