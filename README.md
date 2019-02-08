GENIAL : GENes Identification with Abricate for Lucky biologists
================================================================

Authors : Barbet Pauline, Felten Arnaud

Affiliation: [Food Safety Laboratory - ANSES Maisons Alfort (France)](https://www.anses.fr/en/content/laboratory-food-safety-maisons-alfort-and-boulogne-sur-mer)

You can find the latest version of the tool at [https://github.com/p-barbet/GENIAL](https://github.com/p-barbet/GENIAL)


GENIAL
======

GENIAL aims to identify antimicrobial resistance and virulence genes from bacterial genomes matching them to a database gathering genes of interest using [ABRicate](https://github.com/tseemann/abricate).

### Databases

Default databases available are ([Resfinder](https://cge.cbs.dtu.dk/services/ResFinder/), [CARD](https://card.mcmaster.ca/), [ARG-ANNOT](http://backup.mediterranee-infection.com/article.php?laref=282&titre=arg-annot), [NCBI](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA313047), [EcOH](https://github.com/katholt/srst2/tree/master/data), [PlasmidFinder](https://cge.cbs.dtu.dk/services/PlasmidFinder/), [Ecoli_VF](https://github.com/phac-nml/ecoli_vf) and [VFDB](http://www.mgc.ac.cn/VFs/))

As well as this databes, it's posible to use your own database.

The tool is divided into two scripts.

### GENIALanalysis

GENIALanalysis aims to run ABricate. It takes in input a tsv file containing genomes fasta files paths and IDs.If you want to use your own database you also need to provide a multifasta whith genes IDs as headers. Then the script run ABricate and produce in output one ABRicate result file per genome, corresponding to a tsv file including genes found in each sample.

### GENIALresults

GENIALresults aims to conditionning ABRicate results in the form of matrixes and heatmaps of presence/absence. It takes in input a temporary file produced by the Abricate analysis containing the genomes Abricate results paths and IDs. In the case of vfdb database a file containing the virulence factors names, their family and species is automticaly included in the script.

The output depending on the database used :

* In any cases a matrix in tsv format and a heatmap in png format with all genes found are created


On top of that:

* If you use one of the default databases [Resfinder](https://cge.cbs.dtu.dk/services/ResFinder/) or [VFDB](http://www.mgc.ac.cn/VFs/) news matrix and heatmap by gene type are produced with a correspondace table between the gene name, its family and its number in all genomes.

* If you don't use one of the two previous databases or if you use your own database, only a corespondance table between the gene name and its number in all genomes is produced in addition.

![](workflow.PNG?raw=true "script workflow")


Dependencies
============

The script has been developed with python 3.6 (tested with 3.6.6)

### External dependencies

* [ABRicate](https://github.com/tseemann/abricate) tested with 0.8.7
* [Pandas](https://pandas.pydata.org/) tested with 0.23.4
* [seaborn](https://seaborn.pydata.org/installing.html) tested with 0.9.0


Parameters
==========

### Command line options


|   Options   |                                                              Description                                                              |        Required        |      Default     |
|:-----------:|:-------------------------------------------------------------------------------------------------------------------------------------:|:----------------------:|:----------------:|
|      -f     |                                            tsv file with FASTA files paths ans strains IDs                                            |           Yes          |                  |
|     -dbp    |                                   Path to ABRicate databases repertory. Implies -dbf and --privatedb                                  |   Yes if --privatedb   |                  |
|     -dbf    |                           Multifasta containing the private database sequences. Implies -dbp and --privatedb                          |   Yes if --privatedb   |                  |
|      -T     |                                                        Number of thread to use                                                        |           No           |         1        |
|      -w     |                                                           Working directory                                                           |           No           |         .        |
|      -r     |                                                         Results directory name                                                        |           No           | ABRicate_results |
| --defaultdb |  default databases available : resfinder, card, argannot, acoh, ecoli_vf, plasmidfinder, vfdb or ncbi. Incompatible with --privatedb  | Yes if not --privatedb |                  |
| --privatedb |                              Private database name. Implies -dbp and -dbf. Incompatible with --defaultdb                              | Yes if not --defaultdb |                  |
|   --mincov  |                                                   Minimum proportion of gene covered                                                  |           No           |        80        |
|   --minid   |                                             Minimum proportion of exact nucleotide matches                                            |           No           |        90        |
|     --R     |                                          Remove genes present in all genomes from the matrix                                          |           No           |       False      |

Test 
====

After installing ABRicate and Pandas and seaborn you can test the script with the command line :

## Default database

	python AntiViruce.py -f input_file.tsv --defaultdb vfdb -r results_directory --minid 90 --mincov 80

## Private database

	python AntiViruce.py -f input_file.tsv  --privatedb private_db_name -T 10 -r results_directory --minid 90 --mincov 80 -dbp path_to_abricate_databases_repertory -dbf private_db_multifasta_path







