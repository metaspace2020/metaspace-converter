from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()
long_description = (here / 'README.md').read_text(encoding='utf-8')

setup(
    name='metaspace2anndata',
    version='0.1.0',

    description='Convert metaspace datasets to AnnData',
    long_description=long_description,
    long_description_content_type='text/markdown',

    url='https://git.embl.de/grp-alexandrov/metaspace-to-anndata',
    project_urls={  # Optional
        'Source': 'https://git.embl.de/grp-alexandrov/metaspace-to-anndata',
        #'Publication': ""
    },

    author="Tim Daniel Rose",
    author_email="tim.rose@embl.de",

    license='GPLv3',

    packages=find_packages(),
    install_requires=[
        'pandas',
        'numpy',
        'anndata',
        'metaspace2020'
    ],
    python_requires=">=3.8",

    zip_safe=False,


    classifiers=[
        "License :: OSI Approved :: GNU General Public License v3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Operating System :: MacOS",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3"
    ]
)
