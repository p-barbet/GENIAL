#!/usr/bin/python3
# -*- coding: iso-8859-1 -*-
import os, sys, time
import argparse
import pandas as pd


def get_parser() :
	# Fonction permettant de pourvoir demander des arguments

	parser = argparse.ArgumentParser(description= \
		'Find ARM, virulence and toxin genes running ABRicate')

	parser.add_argument('-f', action="store", dest='input', 
					type=str, required=True, help='tsv file with FASTA files paths and strains IDs (REQUIRED)')

	parser.add_argument('-db', action="store", dest='db_name', \
						type=str, default='resfinder', help='database to use (resfinder, card, \
							argannot, ecoh, ecoli_vf, plasmidfinder, vfdb, ncbi) (default:resfinder)')

	parser.add_argument('-T', action="store", dest='nb_tread', 
					type=str, default=1, help='number of theard to use (default:1)')

	parser.add_argument('-w', action="store", dest='workdir', 
		type=str, default='.', help='working directory (default:current directory)')

	parser.add_argument('-r', action="store", dest='dirres', 
					type=str, default='ABRicate_results', help='results directory name (default:ABRicate_results)')

	parser.add_argument('--mincov', action="store", dest='mincov', \
						type=str, default='80', help='Minimum proportion of gene covered')

	parser.add_argument('--minid', action="store", dest='minid', \
						type=str, default='90', help='Minimum proportion of exact nucleotide matches')

	parser.add_argument('-o', action="store", dest='output_file', \
						type=str, default='ABRicate_files.tsv', help='output file name (default:ABRicate_files.tsv)')

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


def run_ABRicate(fasta_files, db_name, mincov, minid, dir_analysis, nb_threads) :
# Fonction qui lance ABRicate pour chaque souche (par défaut mincov = 80 et minid = 90)

	db_name = db_name.lower() # met le nom de la base de données en minuscules

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

#def fonction_qui_modifie_les_fichiers_abricate() : 


# def get_summaries(abricate_files, dir_summaries) :
# # Fonction qui réalise le sommaire des fichiers résultats d'ABRicate et supprime les fichiers vide

# 	summaries = []

# 	for abricate_file in abricate_files :
# 		summary = dir_summaries + "ABRicate_summary_" + abricate_file + ".tsv" # nom du fichier summary 
# 		os.system("abricate --summary " + abricate_files[abricate_file] + " > " + summary)
# 		summaries.append(summary)


# 		df_summary = pd.read_csv(summary, sep='\t', index_col=0, dtype = str)

# 		files_name = [abricate_file]*len(df_summary.index)
# 		df_summary.index = files_name
# 		df_summary.index.name = "#FILE"

# 		df_summary.to_csv(summary, sep='\t')


# 	return summaries


def get_abricate_files_list(res_dir, abricate_files, output_file) :
# Fonction qui crée un fichier contenant les chemin des fichiers résultats et les IDs des souches

	abricate_list = open(res_dir + output_file, 'w') # ouverture du fichier en écriture

	for abricate_file in abricate_files :
		abricate_list.write(abricate_files[abricate_file] + '\t' + abricate_file + '\n') # écriture du fichier ligne par ligne

	abricate_list.close() # cloture du fichier

	
def main():
	
	##################### gets arguments #####################

	parser=get_parser()
	
	#print parser.help if no arguments
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	
	# mettre tout les arguments dans la variable Argument
	Arguments=parser.parse_args()

	t1 = time.time()

	WORKDIR = Arguments.workdir
	DIRRES = Arguments.dirres
	if WORKDIR[-1] != '/' :
		WORKDIR += '/'

	if DIRRES[-1] != '/' : 
		DIRRES += '/'


	if not os.path.exists(WORKDIR + DIRRES) :
		os.system('mkdir ' + WORKDIR + DIRRES) 
		os.system('mkdir ' + WORKDIR + DIRRES + 'analysis/')

	DIR = WORKDIR + DIRRES
	DIR_ANALYSIS = WORKDIR + DIRRES + "analysis/" # chemin du répertoirequi contiendra les fichier résultats
	
	fasta_files = get_fasta_paths(Arguments.input) # récupère la liste des fichers fasta

	abricate_files = run_ABRicate(fasta_files, Arguments.db_name, Arguments.mincov, Arguments.minid, DIR_ANALYSIS, Arguments.nb_tread) # lance Abricate pour toutes les souches
	#summaries = get_summaries(abricate_files, DIR_SUMMARIES) # Réalise les summary 
	
	get_abricate_files_list(DIR, abricate_files, Arguments.output_file) # créer un fichier avec les liste des fichiers résultats de l'analyse ABricate et les IDs des souches 
	


	t2 = time.time()

	diff = round(t2 - t1,3)
	print ("Temps : " + str(diff) + " secondes")

# lancer la fonction main()  au lancement du script
if __name__ == "__main__":
	main()	            		           		
