#!/usr/bin/env python2

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="GENIAL",
    version="1.0.7",
    author="Pauline Barbet, Arnaud Felten",
    author_email="pauline.barbet@anses.fr, arnaud.felten@anses.fr",
    description="GENIAL: GENes Indentification with Abricate for Lucky biologists",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/p-barbet/GENIAL",
    packages=setuptools.find_packages(),
    classifiers=[
        "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
        "Programming Language :: Python :: 3",
        "Operating System :: POSIX :: Linux",
    ],
        scripts=['src/GENIAL',
             "src/GENIALanalysis",
             "src/GENIALresults",
             "src/GENIALmultidb",
             "src/GENIALupdatedbs",
             "src/GENIALsetupdbs",
             "src/GENIAL_slurm.sh",
             "src/GENIALslurm",
             ],
	data_files = [('GENIALfiles/dbs/enterotox_staph',['dbs/enterotox_staph/sequences']), ('GENIALfiles/dbs/resfinder',['dbs/resfinder/sequences']), ('GENIALfiles/dbs/vfdb',['dbs/vfdb/sequences']), ('GENIALfiles/dbs/phages',['dbs/phages/sequences']), ('GENIALfiles/',['dbs/vfdb/VFs.tsv'])],
    include_package_data=True,
    install_requires=['pandas',   
                      'seaborn',
                      'biopython',
                      'xlrd',
                      ], 
    dependency_links=['https://github.com/tseemann/abricate'], 
    zip_safe=False,



)
