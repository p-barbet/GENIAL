GENIAL : GENes Identification with Abricate for Lucky biologists
================================================================

Authors : Barbet Pauline, Felten Arnaud

Affiliation: [Food Safety Laboratory - ANSES Maisons Alfort (France)](https://www.anses.fr/en/content/laboratory-food-safety-maisons-alfort-and-boulogne-sur-mer)

You can find the latest version of the tool at [https://github.com/p-barbet/GENIAL](https://github.com/p-barbet/GENIAL)


GENIAL
======

GENIAL aims to identify antimicrobial resistance and virulence genes from bacterial genomes matching them to a database gathering genes of interest using [ABRicate](https://github.com/tseemann/abricate).

### Databases

Four default databases are available : ([resfinder](https://cge.cbs.dtu.dk/services/ResFinder/), [vfdb](http://www.mgc.ac.cn/VFs/)), spi and enterotox_staph.

The table below show which kind of research allow each database :


|     Database    |           Research type          |
|:---------------:|:--------------------------------:|
|    resfinder    |  antimicrobial resistance genes  |
|       vfdb      |   virulences and toxins genes    |
|      spi        | Salmonella Pathogenicity Island  |
| enterotox_staph |   staphylococcus enterotoxins    |


As well as this databases, it's posible to use your own database providing in input the sequences in fasta format with gene IDs as headers.

GENIAL is composed by several scripts :

### GENIALanalysis

GENIALanalysis aims to run ABricate. It takes in input a tsv file containing genomes fasta files paths and IDs. If you want to use your own database you also need to provide a multifasta whith genes IDs as headers. Then the script run ABricate and produce in output one ABRicate result file per genome, corresponding to a tsv file including genes found in each sample.

## Command for a default database

	GENIALanalysis -f input_file.tsv -T nb_threads -defaultdb db_name -r res_dir_name -minid min_id -mincov min_cov -w work_dir

## Command for your own database
	
	GENIALanalysis -f input_file.tsv -T nb_threads -privatedb db_name -dbf db_fasta.fasta -r res_dir_name -minid min_id -mincov min_cov -w work_dir


### GENIALresults

GENIALresults aims to conditionning ABRicate results in the form of presence/absence matricies and heatmaps. It takes in input a file produced by Abricateanalysis containing ABRicate results paths and IDs. In the case of vfdb database a file containing the virulence factors names, their family and species is automticaly included in the script.

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

## Command for a default database

	GENIALresults -f ABRicate_files.tsv -w work_dir -r res_dir_name -defaultdb db_name

## Command for your own database

	GENIALresults -f ABRicate_files.tsv -w work_dir -r res_dir_name -privatedb


### GENIAL

GENIAL run the two previous scripts following this workflow for one database :

![](workflow.png?raw=true "script workflow")

## Command for a default database

	GENIAL -f input_file.tsv -T nb_threads -defaultdb db_name -r res_dir_name -minid min_id -mincov min_cov -w work_dir

## Command for your own database

	GENIAL -f input_file.tsv -T nb_threads -w work_dir -r results_directory -minid min_id -mincov min_cov -privatedb db_name -dbf db_fasta.fasta


### GENIALmultidb

GENIALmultidb run GENIAL for several default databases at once and produce a matrix merging all matricies. This script does not take charge of personal databases.

## Command

	GENIALmultidb -f input_file.tsv -db dbs_names -T nb_threads -r res_dir_name -w work_dir -minid min_id -mincov min_cov 


### GENIALupdatedbs

GENIALupdate_databases update resfinder or vfdb databases and setup them in ABRicate.

## Command

	GENIALupdatedbs -db db_name -dbp abricate_db_path


Parameters
==========

### Command line options for GENIAL


|  Options   |                                                              Description                                                              |        Required        |      Default     |
|:----------:|:-------------------------------------------------------------------------------------------------------------------------------------:|:----------------------:|:----------------:|
|     -f     |                                            tsv file with FASTA files paths ans strains IDs                                            |           Yes          |                  |
| -defaultdb |             default databases available : resfinder, vfdb, phages and enterotox staph. Incompatible with --privatedb                  | Yes if not --privatedb |                  |
| -privatedb |                              Private database name. Implies -dbp and -dbf. Incompatible with --defaultdb                              | Yes if not --defaultdb |                  |
|    -dbf    |                           Multifasta containing the private database sequences. Implies -dbp and --privatedb                          |   Yes if --privatedb   |                  |
|     -T     |                                                        Number of thread to use                                                        |           No           |         1        |
|   -mincov  |                                                   Minimum proportion of gene covered                                                  |           No           |        80        |
|   -minid   |                                             Minimum proportion of exact nucleotide matches                                            |           No           |        90        |
|     -w     |                                                           Working directory                                                           |           No           |         .        |
|     -r     |                                                         Results directory name                                                        |           No           | ABRicate_results |


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

You can easily install GENIAL with conda running this commands :

	conda create -n genial
	conda activate genial
	conda config --add channels pbarbet
	conda install genial
	GENIALsetupdbs

After this steps you can use each script described above.

