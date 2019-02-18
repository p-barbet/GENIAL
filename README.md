GENIAL : GENes Identification with Abricate for Lucky biologists
================================================================

Authors : Barbet Pauline, Felten Arnaud

Affiliation: [Food Safety Laboratory - ANSES Maisons Alfort (France)](https://www.anses.fr/en/content/laboratory-food-safety-maisons-alfort-and-boulogne-sur-mer)

You can find the latest version of the tool at [https://github.com/p-barbet/GENIAL](https://github.com/p-barbet/GENIAL)


GENIAL
======

GENIAL aims to identify antimicrobial resistance and virulence genes from bacterial genomes matching them to a database gathering genes of interest using [ABRicate](https://github.com/tseemann/abricate).

### Databases

Four default databases are available : [resfinder](https://cge.cbs.dtu.dk/services/ResFinder/), [vfdb](http://www.mgc.ac.cn/VFs/)), phages(https://www.ebi.ac.uk/genomes/phage.html) and enterotox_staph.

The table below show which kind of research allow each database :


|     Database    |          Research type         |
|:---------------:|:------------------------------:|
|    resfinder    | antimicrobial resistance genes |
|       vfdb      |   virulences and toxins genes  |
|      phages     |             phages             |
| enterotox_staph |   staphylococcus enterotoxins  |


As well as this databases, it's posible to use your own database providing in input the sequences in fasta format with gene IDs as headers.

GENIAL is composed of several scripts :

### GENIALanalysis

GENIALanalysis aims to run ABricate. It takes in input a tsv file containing genomes fasta files paths and IDs. If you want to use your own database you also need to provide a multifasta whith genes IDs as headers. Then the script run ABricate and produce in output one ABRicate result file per genome, corresponding to a tsv file including genes found in each sample.


### GENIALresults

GENIALresults aims to conditionning ABRicate results in the form of presence/absence matricies and heatmaps. It takes in input a temporary file produced by the Abricate analysis containing the genomes Abricate results paths and IDs. In the case of vfdb database a file containing the virulence factors names, their family and species is automticaly included in the script.

The output depending on the database used as presented below :

| Legend |                    Possible outputs                    |
|:------:|:------------------------------------------------------:|
|    1   |               Matrix and heatmap by genes              |
|    2   |          Matrix and heatmap by genes families          |
|    3   |       Correspondence table gene name/gene number       |
|    4   | Correspondence table gene name/gene family/gene number |


|     Database     |  Outputs  |
|:----------------:|:---------:|
|       vfdb       | 1 + 2 + 3 |
|     resfinder    | 1 + 2 + 3 |
|      phages      |   1 + 4   |
|  enterotox_staph |   1 + 4   |
| private database |   1 + 4   |


Matricies are produced in tsv format and a heatmaps in png format.


### GENIAL

GENIAL run the two previous scripts following this workflow for one database :

![](workflow.PNG?raw=true "script workflow")


### GENIALmultidb

GENIALmultidb run GENIAL for several databases and produce a matrix merging all matricies.


### GENIALupdatedbs

GENIALupdate_databases update resfinder or vfdb databases in setup them in ABRicate.


Parameters
==========

### Command line options for GENIAL


|   Options   |                                                              Description                                                              |        Required        |      Default     |
|:-----------:|:-------------------------------------------------------------------------------------------------------------------------------------:|:----------------------:|:----------------:|
|      -f     |                                            tsv file with FASTA files paths ans strains IDs                                            |           Yes          |                  |
|     -dbp    |                                   Path to ABRicate databases repertory. Implies -dbf and --privatedb                                  |   Yes if --privatedb   |                  |
|     -dbf    |                           Multifasta containing the private database sequences. Implies -dbp and --privatedb                          |   Yes if --privatedb   |                  |
|      -T     |                                                        Number of thread to use                                                        |           No           |         1        |
|      -w     |                                                           Working directory                                                           |           No           |         .        |
|      -r     |                                                         Results directory name                                                        |           No           | ABRicate_results |
| --defaultdb |             default databases available : resfinder, vfdb, phages and enterotox staph. Incompatible with --privatedb                  | Yes if not --privatedb |                  |
| --privatedb |                              Private database name. Implies -dbp and -dbf. Incompatible with --defaultdb                              | Yes if not --defaultdb |                  |
|   --mincov  |                                                   Minimum proportion of gene covered                                                  |           No           |        80        |
|   --minid   |                                             Minimum proportion of exact nucleotide matches                                            |           No           |        90        |
|     --R     |                                          Remove genes present in all genomes from the matrix                                          |           No           |       False      |


Dependencies
============

GENIAL has been developed with python 3.6 (tested with 3.6.6)

### External dependencies

* [ABRicate](https://github.com/tseemann/abricate) tested with 0.8.7
* [pandas](https://pandas.pydata.org/) tested with 0.23.4
* [seaborn](https://seaborn.pydata.org/installing.html) tested with 0.9.0
* [biopython](https://biopython.org/wiki/Download) tested with 1.72


Installation
============

To install GENIAL run this command lines :

	conda config --add channels pbarbet
	conda install genial

After installing GENIAL you can use simply each script typing their name in a terminal.


Test 
====

## Default database

	python GENIAL -f input_file.tsv --defaultdb vfdb -r results_directory --minid 90 --mincov 80

## Private database

	python GENIAL -f input_file.tsv  --privatedb private_db_name -T 10 -r results_directory --minid 90 --mincov 80 -dbp path_to_abricate_databases_repertory -dbf private_db_multifasta_path