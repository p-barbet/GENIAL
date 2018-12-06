#!/usr/bin/python3
# -*- coding: iso-8859-1 -*-
import os, sys, time
import argparse
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from pandas.tools.plotting import table
import imgkit
import copy


def get_parser() :
	# Fonction permettant de pourvoir demander des arguments

	parser = argparse.ArgumentParser(description= \
		'Make a presence/absence matrix with ABRicate results files')

	parser.add_argument('-f', action="store", dest='abricate_list', 
					type=str, required=True, help='summaries files from ABRicate (path)(REQUIRED)')

	parser.add_argument('-w', action="store", dest='workdir', 
					type=str, default='.', help='working directory (default:current directory)')

	parser.add_argument('-r', action="store", dest='dirres', 
					type=str, default='ABRicate_results', help='results directory name (default:ABRicate_results)')

	parser.add_argument('--research', dest='research', choices=['antibiotic_resistance','virulence_factor'], required=True,
					help='Type of researched genes antibiotic antibiotic resistance or virulence factor (')

	parser.add_argument('--R', action="store_true", dest='remove',
					default=False, help='remove genes present in all genomes from the matrix (default:False)')

	return parser


def get_paths_files(input_file) :
# Fonction qui récupère la liste des fichiers résultats d'ABRicate et les IDs des souches
	
	dico_abricate = {}

	data = open(input_file, 'r')
	lines = data.readlines() 
	data.close()

	for line in lines :

		line = line.rstrip() # retire les retours chariot
		infos = line.split("\t") # split chaque ligne selon les tabulations 
		dico_abricate[infos[1]] = infos[0] # remplis le dico avec les IDs en clé et la chemins vers les fichiers résultats d'ABRicate en valeurs

	return dico_abricate


def get_VF_family_dico(family_file) :

	data = open(family_file, 'r')
	lines  = data.readlines()
	data.close()

	dico_family = {}

	for line in lines :

		line = line.rstrip()
		infos = line.split('\t')
	

		VF_name = infos[0]
		species_name = infos[1]

		try :
			family = infos[2].split('; ')[0]

		except IndexError:
			family = 'Unclassified'

		try :
			dico_family[species_name][VF_name] = family

		except KeyError:
			dico_family[species_name] = {}
			dico_family[species_name][VF_name] = family

	return dico_family


def get_matrix_with_genes(abricate_files, research_type) :
# Fonction qui réalise la matrice de présence/abcence de tous les gènes pour chaque génome
	
	dfs = [] # liste qui contiendra les datframes crées à partir des fichier ABRicate

	# création d'un dataframe pour chaque fichier ABRicate
	for abricate_file in abricate_files : 

		genes = []
		

		df_abricate = pd.read_csv(abricate_files[abricate_file], sep='\t', index_col=0, dtype = str) # lecture du fichier abricate
		
		genes_list = df_abricate['GENE'] # liste des gènes
		
		ID = [abricate_file] # ID de la souche	

		if len(genes_list) != 0 : # si au moins à gène de résiatnce à été trouvé

			if research_type == "virulence_factor" :
				
				VF_names = []
				species_names = []

				product_infos = df_abricate['PRODUCT']
				
				for info in product_infos :
					VF_name = '|' + info.split('[')[1].split(' (')[0]
					VF_names.append(VF_name)

			
					species_name = '|' + " ".join(info.split('[')[2].split(' ')[0:2])

					if species_name[-1] == "]" :
						species_name = species_name.replace("]", "")

					species_names.append(species_name)
					
				genes_list = genes_list + VF_names + species_names

			for gene in genes_list :
				genes.append(gene) # ajout de chaque gène à la liste gène

			df_genes = pd.DataFrame(1, index = ID, columns = genes) # crée un dataframe d'une ligne remplis de 1 avec les gènes de la souche en colonnes et le nom de la souche en index
			df_genes = df_genes.groupby(df_genes.columns, axis=1).sum()# si un gène à été trouvé pusieurs fois dans le génome fait la somme des colonnes avec le même nom de gene

		else :

			df_genes = pd.DataFrame(index = ID) # sinon création d'un dataframe vide avec comme index l'ID

			
		dfs.append(df_genes) # ajout du dataframe à la liste

	# merge des dataframes
	matrix = pd.concat(dfs, sort=True)

	# remplace les valeurs NaN par des 0
	matrix.fillna(0, inplace=True)

	return matrix


def get_matrix_by_genes_types(MATRIX, research_type) :
# Fonction qui réalise la matrice par famille d'antibiotique

	columns = MATRIX.columns # liste des nom de colonne
	genomes = MATRIX.index # liste des IDs

	if research_type == "virulence_factor" :
		family_names = []
		VF_family_dico = get_VF_family_dico('VFs.csv')

		for column in columns :
			VF_name = column.split('|')[1]
			print(VF_name)

			species = column.split('|')[2]

			try :
				family_name = VF_family_dico[species][VF_name]

			except KeyError :
				family_name = 'Unclassified'

			if not family_name in family_names : 
				family_names.append(family_name)

		matrix_by_VF_family = pd.DataFrame(0, index = genomes, columns = family_names)

		for genome in genomes : # pour chque génome (lignes)
			for column in columns : # pour chaque colonne
				VF_name = column.split('|')[1]
				species = column.split('|')[2]

				try :
					family_name = VF_family_dico[species][VF_name]

				except KeyError :
					family_name = 'Unclassified'

				if MATRIX[column][genome] != 0 : # si le gène est présent dans le génome
					matrix_by_VF_family[family_name][genome] += MATRIX[column][genome]  # incrémente la matrice de sa valeur

		return matrix_by_VF_family

	else :
		antibiotics = [] # liste des familles d'antibio
	
		for column in columns :
			antibiotic = column.split("|")[-1] # nom de la famille de l'antibiotique à laquelle le gène est résistant

			if not antibiotic in antibiotics : 
					antibiotics.append(antibiotic) # si la famille n'est pas déjà présente dans la liste, elle  y est ajoutée
	
		matrix_by_antibiotics = pd.DataFrame(0, index = genomes, columns = antibiotics) # nouvelle matrice remplie de 0 avec les génomes en index et les familles d'antibiotique en colonne

		for genome in genomes : # pour chque génome (lignes)
			for column in columns : # pour chaque colonne
				antibiotic = column.split("|")[-1] # nom de la famille d'antibiotique

				if MATRIX[column][genome] != 0 : # si le gène est présent dans le génome
					matrix_by_antibiotics[antibiotic][genome] += MATRIX[column][genome]  # incrémente la matrice de sa valeur

		return matrix_by_antibiotics


def write_matrix(MATRIX, matrix_name, dir_matrix, index) :
# Fonction qui écrit la matrice dans un fichier

	MATRIX.to_csv(dir_matrix + matrix_name, sep='\t', index = index) # écriture de la matrice dans un fichier tsv


def get_correspondance_table(MATRIX, research_type) :
# Fonction qui réalise une table de correspondance en attribuant à chaque gène son type de résistance et le nombre de fois qu'il ets retrouvé tous les génomes confondu

	columns = MATRIX.columns # nom des colonnes

	if research_type == "virulence_factor" :
		genes = []
		VF_names = []
		family_names = []
		species_names = []
		number = []

		VF_family_dico = get_VF_family_dico('VFs.csv')

		cor_table = pd.DataFrame(columns = ['genes', 'species', 'VF names', 'family names', 'number']) # Définition d'un dataframe vide avec 3 colonnes

		for column in columns : # pour chaque colonne de la matrice

			gene = column.split("|")[0] # nom du gène

			if not gene in genes : # si le gène n'est pas dans la liste
				genes.append(gene) # ajout du gène à la liste
				number.append(sum(MATRIX[column])) # exemplaire du gène dans tous les génomes

				VF_name = column.split('|')[1]
				VF_names.append(VF_name)

				species_name = column.split('|')[2]
				species_names.append(species_name)

				try :
					family_name = VF_family_dico[species_name][VF_name]

				except KeyError :
					family_name = 'Unclassified'

				family_names.append(family_name)

		cor_table['genes'] = genes # remplissage de la colonne gene
		cor_table['species'] = species_names
		cor_table['VF names'] = VF_names # remplissage de la colonne antibiotique
		cor_table['family names'] = family_names
		cor_table['number'] = number # remplissage de la colonne effectif

		# if all(MATRIX[VF names] == ) : # Si toutes les valuers de la colonne du gène sout égales à 1
		# 	MATRIX.pop(gene) # suppression du gène


	else :

		genes = [] # liste des gènes
		antibiotics = [] # liste des antibiotiques
		number = [] # exemplaire du gène dans tous les génomes

		cor_table = pd.DataFrame(columns = ['genes', 'antibiotic resistance', 'number']) # Définition d'un dataframe vide avec 3 colonnes
	
		for column in columns : # pour chaque colonne de la matrice
			gene = column.split("|")[0] # nom du gène

			if not gene in genes : # si le gène n'est pas dans la liste
				genes.append(gene) # ajout du gène à la liste

				number.append(sum(MATRIX[column])) # exemplaire du gène dans tous les génomes

				antibiotic = column.split("|")[-1] # antibiotique
				antibiotics.append(antibiotic) # 

		cor_table['genes'] = genes # remplissage de la colonne gene
		cor_table['antibiotic resistance'] = antibiotics # remplissage de la colonne antibiotique
		cor_table['number'] = number # remplissage de la colonne éffectif

	print(cor_table)
	print(cor_table['number'].sum())

	return cor_table


def remove_genes_types_from_genes(MATRIX) : 
# Fonction qui renomme les collones de la matrice avec tous les gène identifiés

	columns = MATRIX.columns # nom des colonnes

	genes = [] # liste des gènes

	for column in columns : # Pour chaque clonnes
		gene = column.split('|')[0] # nom du gène
		genes.append(gene) # ajout à la liste des gènes

	MATRIX.columns = genes # attribution de la liste aux colonnes

	return MATRIX


def remove_genes_present_in_all_genomes(MATRIX) :
# Fonction qui supprime dans la matrice les gènes présent dans tous les génomes et en fait le liste
	genes = MATRIX.columns # nom des gènes

	removes_genes = [] # liste des gènes supprimés

	for gene in genes : # pour chaque gène
		print(gene)
		print(all(MATRIX[gene] == 1)) 
		if all(MATRIX[gene] == 1) : # Si toutes les valuers de la colonne du gène sout égales à 1
			MATRIX.pop(gene) # suppression du gène
			removes_genes.append(gene) # ajout du gène à la liste des gènes supprimés

	print(removes_genes)

	return MATRIX


def heatmap(MATRIX, heatmap_name, dir_heatmap) :
# fonction qui réalise une heatmap

	# sns.clustermap(df, metric="correlation", method="single", cmap="Blues", standard_scale=1, row_colors=row_colors)
	
	try : # si pas de RecursionError

		heatmap = sns.clustermap(MATRIX, metric="euclidean", method="average", figsize=(18, 14), linewidths=.003) # réalisation de la heatmap
		heatmap.savefig(dir_heatmap + heatmap_name) # sauvegarde de la heatmap dans un png

	except RecursionError: # sinon

		print("Erreur lors de la construction de la heatmap, la matrice est trop grande") # affichage d'un massage d'erreur

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
		os.system('mkdir ' + WORKDIR + DIRRES + 'matrix/')
		os.system('mkdir ' + WORKDIR + DIRRES + 'heatmap/')

	else : 
		os.system('mkdir ' + WORKDIR + DIRRES + 'matrix/')
		os.system('mkdir ' + WORKDIR + DIRRES + 'heatmap/')


	DIR_MATRIX = WORKDIR + DIRRES + "matrix/" # chemin de destination des matrices
	DIR_HEATMAP = WORKDIR + DIRRES + "heatmap/" #chemin de derstinatin des heatmap


	abricate_files = get_paths_files(Arguments.abricate_list) # récupère les fichier résulatant de l'analyse ABRicate
	
	matrix_all_genes = get_matrix_with_genes(abricate_files, Arguments.research) # réalise la matrice

	matrix_by_genes_types = get_matrix_by_genes_types(matrix_all_genes, Arguments.research)

	cor_table = get_correspondance_table(matrix_all_genes, Arguments.research)

	matrix_all_genes = remove_genes_types_from_genes(matrix_all_genes)

	if Arguments.remove : 

		matrix_all_genes = remove_genes_in_all_genomes(matrix_all_genes)

	write_matrix(matrix_all_genes, "matrix_all_genes.tsv", DIR_MATRIX, True) # ecrit la matrice des gènes
	write_matrix(matrix_by_genes_types, "matrix__by_gene_type.tsv", DIR_MATRIX, True) # ecrit la matrice des antibio
	write_matrix(cor_table, "correspondance_table.tsv", DIR_MATRIX, False) # ecrit la table de correspondance
	


	heatmap(matrix_all_genes, "heatmap_all_genes.png", DIR_HEATMAP) # réalise la heatmap des gènes
	heatmap(matrix_by_genes_types, "heatmap_by_gene_type.png", DIR_HEATMAP) # réalise la heatmap des antibio

	


	t2 = time.time()

	diff = round(t2 - t1,3)
	print ("Temps : " + str(diff) + " secondes")

# lancer la fonction main()  au lancement du script
if __name__ == "__main__":
	main()	  
