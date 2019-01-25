#!/usr/bin/python3
# -*- coding: utf-8 -*-

import os, sys, time
import argparse

def get_parser() :
	# Fonction permettant de pourvoir demander des arguments

	parser = argparse.ArgumentParser(description= \
		'Find ARM, virulence and toxin genes running ABRicate', \
			formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	parser.add_argument('-f', action="store", dest='input', 
					type=str, required=True, help='tsv file with FASTA files paths and strains IDs (REQUIRED)')

	db_type = parser.add_mutually_exclusive_group(required=True)

	db_type.add_argument('--defaultdb', action="store", dest='default_database', \
						type=str, choices=['resfinder', 'card',	'argannot', 'ecoh', \
							'ecoli_vf', 'plasmidfinder', 'vfdb', 'ncbi', 'vir_clost', 'enterotox_staph'], help='default \
								database to use (resfinder, card, argannot, ecoh, ecoli_vf, plasmidfinder, vfdb, ncbi. Incompatible with --privatedb)')

	db_type.add_argument('--privatedb', action="store", dest='private_database', \
						type=str, help='private database name. Implies -dbp, -dbf. Incompatible with --defaultdb')

	parser.add_argument('-dbp', action="store", dest='private_db_path', \
						type=str, help='path to abricate \
							databases repertory. Implies -dbf, --privatedb')	

	parser.add_argument('-dbf', action="store", dest='private_db_fasta', \
						type=str, help='Multifasta containing \
							the private database sequences. Implies -dbp, --privatedb')	

	parser.add_argument('-T', action="store", dest='nb_thread', 
					type=str, default=1, help='number of theard to use')

	parser.add_argument('-w', action="store", dest='workdir', 
		type=str, default='.', help='working directory')

	parser.add_argument('-r', action="store", dest='resdir', 
					type=str, default='ABRicate_results', help='results directory name')

	parser.add_argument('--mincov', action="store", dest='mincov', \
						type=str, default='80', help='Minimum proportion of gene covered')

	parser.add_argument('--minid', action="store", dest='minid', \
						type=str, default='90', help='Minimum proportion of exact nucleotide matches')

	parser.add_argument('--R', action="store_true", dest='remove',
					default=False, help='remove genes present in all genomes from the matrix (default:False)')

	return parser








	

def main():
	
	##################### gets arguments #####################

	parser=get_parser()
	
	#print parser.help if no arguments
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	
	# mettre tout les arguments dans la variable Argument
	Arguments=parser.parse_args()

	if Arguments.private_database is not None and (Arguments.private_db_path is None or Arguments.private_db_fasta is None) :
		 parser.error("--privatedb argument requires -dbp and -dbf.")

	if Arguments.default_database is not None and (Arguments.private_db_path is not None or Arguments.private_db_fasta is not None) :
		 parser.error("--defaultdb argument not requires -dbp and -dbf.")


	if Arguments.workdir[-1] != '/' :
		Arguments.workdir += '/'

	if Arguments.resdir[-1] != '/' : 
		Arguments.resdir += '/'


	if Arguments.default_database is not None :

		os.system("python abricate_analysis.py -f " +  Arguments.input + " -T " + Arguments.nb_thread + " --defaultdb " + Arguments.default_database + ' -r ' + Arguments.resdir + " --minid " + Arguments.minid + " --mincov " + Arguments.mincov + " -w " + Arguments.workdir)

		if os.path.isfile(Arguments.workdir + Arguments.resdir + "ABRicate_files.tsv") :  # Vérifie que le fichier abricate_files.tsv existe

			if Arguments.remove :
				os.system("python abricate_matrix.py -f " + Arguments.workdir + Arguments.resdir + "ABRicate_files.tsv -w " + Arguments.workdir + " -r " + Arguments.resdir + " --defaultdb " + Arguments.default_database + " --R")

			else :
				os.system("python abricate_matrix.py -f " + Arguments.workdir + Arguments.resdir + "ABRicate_files.tsv -w " + Arguments.workdir + " -r " + Arguments.resdir + " --defaultdb " + Arguments.default_database)


	elif Arguments.private_database is not None :

		os.system("python abricate_analysis.py -f " +  Arguments.input + " -T " + Arguments.nb_thread + " --privatedb " + Arguments.private_database + " -dbp " + Arguments.private_db_path + " -dbf " + Arguments.private_db_fasta + ' -r ' + Arguments.resdir + " --minid " + Arguments.minid + " --mincov " + Arguments.mincov + " -w " + Arguments.workdir)

		if os.path.isfile(Arguments.workdir + Arguments.resdir + "ABRicate_files.tsv") : 

			if Arguments.remove :
				os.system("python abricate_matrix.py -f " + Arguments.workdir + Arguments.resdir + "ABRicate_files.tsv -w " + Arguments.workdir + " -r " + Arguments.resdir + " --privatedb --R")

			else :
				os.system("python abricate_matrix.py -f " + Arguments.workdir + Arguments.resdir + "ABRicate_files.tsv -w " + Arguments.workdir + " -r " + Arguments.resdir + " --privatedb")



	
	

# lancer la fonction main()  au lancement du script
if __name__ == "__main__":
	main()	