#!/usr/bin/python3
# -*- coding: utf-8 -*-
import os, sys, time
import argparse
import pandas as pd


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
							'ecoli_vf', 'plasmidfinder', 'vfdb', 'ncbi'], help='default \
								database to use (resfinder, card, argannot, ecoh, ecoli_vf, plasmidfinder, vfdb, ncbi. Incompatible with --privatedb)')

	db_type.add_argument('--privatedb', action="store", dest='private_database', \
						type=str, help='private database name. Implies -dbp, -dbf. Incompatible with --defaultdb')

	parser.add_argument('-dbp', action="store", dest='private_db_path', \
						type=str, help='path to abricate \
							databases repertory. Implies -dbf, --privatedb')	

	parser.add_argument('-dbf', action="store", dest='private_db_fasta', \
						type=str, help='Multifasta containing \
							the private database sequences. Implies -dbp, --privatedb')	

	parser.add_argument('-T', action="store", dest='nb_tread', 
					type=str, default=1, help='number of theard to use')

	parser.add_argument('-w', action="store", dest='workdir', 
		type=str, default='.', help='working directory')

	parser.add_argument('-r', action="store", dest='resdir', 
					type=str, default='ABRicate_results', help='results directory name')

	parser.add_argument('--mincov', action="store", dest='mincov', \
						type=str, default='80', help='Minimum proportion of gene covered')

	parser.add_argument('--minid', action="store", dest='minid', \
						type=str, default='90', help='Minimum proportion of exact nucleotide matches')

	parser.add_argument('-o', action="store", dest='output_file', \
						type=str, default='ABRicate_files.tsv', help='output file name')

	return parser


def get_fasta_paths(input_file) :
# Fonction qui récupère la liste des fichiers fasta à analyser et l'assoccie aux IDs
	
	dico_fasta = {} 

	# Récupèreation du contenu du fichier
	data = open(input_file, 'r')
	lines = data.readlines() 
	data.close()

	for line in lines :

		line = line.rstrip() # retire les retours chariot
		infos = line.split("\t") # split chaque ligne selon les tabulations 
		dico_fasta[infos[1]] = infos[0] # remplis le dico avec les IDs en clé et la chemins vers les FASTA en valeurs

	return dico_fasta


def setup_private_db(db_name, abricate_dbs_repertory, db_multifasta) :

	if abricate_dbs_repertory[-1] != '/' :
		abricate_dbs_repertory += '/'

	if not os.path.exists(abricate_dbs_repertory + db_name) :
		os.system("mkdir " + abricate_dbs_repertory + db_name)
	else : 
		sys.exit("La base de donnée " + db_name + " existe déjà")

	os.system("cp " + db_multifasta + " " + abricate_dbs_repertory + db_name)

	os.system("abricate --setupdb")


def run_ABRicate(fasta_files, db_name, mincov, minid, dir_analysis, nb_threads) :
# Fonction qui lance ABRicate pour chaque souche (par défaut mincov = 80 et minid = 90)

	abricate_files = {} 

	for fasta_file in fasta_files :

		abricate_result = dir_analysis + "ABRicate_" + db_name + "_" + fasta_file + ".tsv" # nom du fichier contenant le résultat de l'analyse (fasta_file = ID de la souche)
		os.system("abricate " + fasta_files[fasta_file] + " --db " + db_name + " --mincov " + mincov + " --minid " + minid + " --threads " + str(nb_threads) + " > " + abricate_result) # ligne de commade d'ABRicate
		abricate_files[fasta_file] = abricate_result #remplis le dico avec les IDs en clés et la chemins vers les fichier résultats d'ABRicate en valeurs 

		df_result = pd.read_csv(abricate_result, sep='\t', index_col=0, dtype = str) # lecture du fichier résultat avec pandas (dataframe)

		files_name = [fasta_file]*len(df_result.index) # liste contenant autant de fois l'ID de la souche qu'il n'y a de gène dansle fichier
		df_result.index = files_name # renomage du nom du fichier par l'ID de la souche pour tous les gènes
		df_result.index.name = "#FILE" # nom de l'index

		df_result.to_csv(abricate_result, sep='\t') # Récriture du fichier

	return abricate_files


def get_abricate_files_list(res_dir, abricate_files, output_file) :
# Fonction qui crée un fichier contenant les chemin des fichiers résultats et les IDs des souches

	abricate_list = open(res_dir + output_file, 'w') # ouverture du fichier en écriture

	for abricate_file in abricate_files :
		abricate_list.write(abricate_files[abricate_file] + '\t' + abricate_file + '\n') # écriture du fichier ligne par ligne

	abricate_list.close() # cloture du fichier

def uninstall_private_db(db_name, abricate_dbs_repertory) :

	if abricate_dbs_repertory[-1] != '/' :
		abricate_dbs_repertory += '/'
	
	os.system("rm -R " + abricate_dbs_repertory + db_name)

	os.system("abricate --setupdb")

	
def main():
	
	##################### gets arguments #####################

	parser=get_parser()
	
	#print parser.help if no arguments
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	
	# mettre tout les arguments dans la variable Argument
	Arguments=parser.parse_args()

	print(Arguments)

	if Arguments.private_database is not None and (Arguments.private_db_path is None or Arguments.private_db_fasta is None) :
		 parser.error("--privatedb argument requires -dbp and -dbf.")

	print(Arguments.default_database is not None)

	if Arguments.default_database is not None and (Arguments.private_db_path is not None or Arguments.private_db_fasta is not None) :
		 parser.error("--defaultdb argument not requires -dbp and -dbf.")

	t1 = time.time()

	WORKDIR = Arguments.workdir
	RESDIR = Arguments.resdir
	if WORKDIR[-1] != '/' :
		WORKDIR += '/'

	if RESDIR[-1] != '/' : 
		RESDIR += '/'


	if not os.path.exists(WORKDIR + RESDIR) :
		os.system('mkdir ' + WORKDIR + RESDIR) 
		os.system('mkdir ' + WORKDIR + RESDIR + 'analysis/')

	DIR = WORKDIR + RESDIR
	DIR_ANALYSIS = WORKDIR + RESDIR + "analysis/" # chemin du répertoirequi contiendra les fichier résultats
	
	fasta_files = get_fasta_paths(Arguments.input) # récupère la liste des fichers fasta

	setup_private_db(Arguments.private_database, Arguments.private_db_path, Arguments.private_db_fasta)

	#abricate_files = run_ABRicate(fasta_files, Arguments.default_database, Arguments.mincov, Arguments.minid, DIR_ANALYSIS, Arguments.nb_tread) # lance Abricate pour toutes les souches
	
	#get_abricate_files_list(DIR, abricate_files, Arguments.output_file) # créer un fichier avec les liste des fichiers résultats de l'analyse ABricate et les IDs des souches 
	
	uninstall_private_db(Arguments.private_database, Arguments.private_db_path)

	t2 = time.time()

	diff = round(t2 - t1,3)
	print ("Temps : " + str(diff) + " secondes")

# lancer la fonction main()  au lancement du script
if __name__ == "__main__":
	main()	            		           		
