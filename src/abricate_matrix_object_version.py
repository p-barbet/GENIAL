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
		'Make a presence/absence matrix with ABRicate results files',\
			formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	parser.add_argument('-f', action="store", dest='abricateList', 
					type=str, required=True, help='summaries files from ABRicate (path)(REQUIRED)')

	parser.add_argument('-w', action="store", dest='workdir', 
					type=str, default='.', help='working directory (default:current directory)')

	parser.add_argument('-r', action="store", dest='resdir', 
					type=str, default='ABRicate_results', help='results directory name (default:ABRicate_results)')

	db_type = parser.add_mutually_exclusive_group(required=True)

	db_type.add_argument('--defaultdb', action="store", dest='defaultDatabase', \
						type=str, choices=['resfinder', 'card',	'argannot', 'ecoh', \
							'ecoli_vf', 'plasmidfinder', 'vfdb', 'ncbi'], help='default \
								database to use (resfinder, card, argannot, ecoh, ecoli_vf, plasmidfinder, vfdb, ncbi. Incompatible with --privatedb)')

	db_type.add_argument('--privatedb', action="store_true", dest='privateDatabase', \
						help='private database name. Implies -dbp, -dbf. Incompatible with --defaultdb')

	# parser.add_argument('--research', dest='research', choices=['antibiotic_resistance','virulence_factor'], required=True,
	# 				help='Type of researched genes antibiotic antibiotic resistance or virulence factor ')

	parser.add_argument('--R', action="store_true", dest='remove',
					default=False, help='remove genes present in all genomes from the matrix (default:False)')

	return parser

class genome(object) :

	def __init__(self) :
		self.ID = ""
		self.abricateFile = ""
		self.genes = [] #liste des gènes 
		self.abricateMatrix = None # dataframe contenant les informations du fichier
		self.genesMatrix = None


	def setID(self, ID) :
		self.ID = ID

	def setAbricateFile(self, abricateFile) :
		self.abricateFile = abricateFile

	def setAbricateMatrix(self) :
		self.abricateMatrix = pd.read_csv(self.abricateFile, sep='\t', index_col=0, dtype = str) # dataframe contenant les informations du fichie

	def setGenes(self) :

		genesList = self.abricateMatrix['GENE']

		if len(genesList) != 0 : # si au moins à gène de résiatnce à été trouvé
		
			for geneName in genesList :
				self.genes.append(geneName) # ajout de chaque gène à la liste gène

	def setGenesMatrix(self) :
		if len(self.genes) != 0 : # si au moins à gène de résiatnce à été trouvé
		
			for geneName in self.genes :
				self.genesMatrix = pd.DataFrame(1, index = [self.ID], columns = self.genes) # crée un dataframe d'une ligne remplis de 1 avec les gènes de la souche en colonnes et le nom de la souche en index
				self.genesMatrix = self.genesMatrix.groupby(self.genesMatrix.columns, axis=1).sum()# si un gène à été trouvé pusieurs fois dans le génome fait la somme des colonnes avec le même nom de gene

		else :

			self.genesMatrix = pd.DataFrame(index = [self.ID]) # sinon création d'un dataframe vide avec comme index l'ID


	def getGenesObjects(self, dicoGenes, database, dicoVfFamilies) :
		for geneName in self.genes :
			if geneName not in dicoGenes :
				dicoGenes[geneName] = gene()

				productInfos = self.abricateMatrix.set_index('GENE')['PRODUCT'][geneName]

				dicoGenes[geneName].setName(geneName)
				dicoGenes[geneName].setDatabase(database)

				if database == 'vfdb' : 
					dicoGenes[geneName].setVfName(productInfos)
					dicoGenes[geneName].setSpecies(productInfos)
					dicoGenes[geneName].setVfFamily(dicoVfFamilies)

				elif database == 'resfinder' :
					dicoGenes[geneName].setAntibioticFamily(productInfos)

	

class gene(object) :

	def __init__(self) :
		self.name = ""
		self.database = ""

		# si vfdb
		self.vfName = ""
		self.vfFamily = ""
		self.species = ""

		# si resfinder
		self.antibioticFamily = ""

	def setName(self, geneName) :
		self.name = geneName

	def setDatabase(self, database) :
		self.database = database

	def setVfName(self, productInfos) :
		self.vfName = productInfos.split('[')[1].split(' (')[0]
			
	def setSpecies(self, productInfos) :
		species = " ".join(productInfos.split('[')[2].split(' ')[0:2])

		if species[-1] == "]" :
			species = species.replace("]", "")

		self.species = species

	def setVfFamily(self, dicoVfFamilies) :
		if (self.species in dicoVfFamilies) and (self.vfName in dicoVfFamilies[self.species]) :
			self.vfFamily = dicoVfFamilies[self.species][self.vfName]

		else : 
			self.vfFamily = 'Unclassified'

	def setAntibioticFamily(self, productInfos) :
		self.antibioticFamily = productInfos




def getGenomesObjects(inputFile, dicoGenomes, database, dicoGenes, dicoVfFamilies) :
	
	data = open(inputFile, 'r')
	lines = data.readlines() 
	data.close()

	for line in lines :

		line = line.rstrip() # retire les retours chariot
		infos = line.split("\t") # split chaque ligne selon les tabulations 

		ID = infos[1]
		abricateFile = infos[0]

		dicoGenomes[ID] = genome()

		dicoGenomes[ID].setID(ID)
		dicoGenomes[ID].setAbricateFile(abricateFile)
		dicoGenomes[ID].setAbricateMatrix()
		dicoGenomes[ID].setGenes()
		dicoGenomes[ID].setGenesMatrix()


		dicoGenomes[ID].getGenesObjects(dicoGenes, database, dicoVfFamilies)



def getVfFamiliesDico(familiesFile) :

	data = open(familiesFile, 'r')
	lines  = data.readlines()
	data.close()

	dicoVfFamilies = {}

	for line in lines :

		line = line.rstrip()
		infos = line.split('\t')

		vfName = infos[0]
		species = infos[1]

		try :
			vfFamily = infos[2].split('; ')[0]

		except IndexError:
			vfFamily = 'Unclassified'

		if species in dicoVfFamilies :
			dicoVfFamilies[species][vfName] = vfFamily

		else :
			dicoVfFamilies[species] = {}
			dicoVfFamilies[species][vfName] = vfFamily

	return dicoVfFamilies



def getMatrixAllGenes(dicoGenomes) :
	dataframes = []

	for genome in dicoGenomes :
		dataframes.append(dicoGenomes[genome].genesMatrix)

	# merge des dataframes
	matrixAllGenes = pd.concat(dataframes, sort=True)

	# remplace les valeurs NaN par des 0
	matrixAllGenes.fillna(0, inplace=True)

	return matrixAllGenes



def getVfdbMatrixByGenesTypes(matrixAllGenes, dicoGenes) :
	genomes = matrixAllGenes.index
	genesList = matrixAllGenes.columns

	vfFamilies = []

	for geneName in genesList :
		vfFamily = dicoGenes[geneName].vfFamily

		if vfFamily not in vfFamilies :
			vfFamilies.append(vfFamily)

	matrixByVfFamilies = pd.DataFrame(0, index = genomes, columns = vfFamilies)

	for genome in genomes :
		for geneName in genesList :
			familyName = dicoGenes[geneName].vfFamily

			if matrixAllGenes[geneName][genome] != 0 : # si le gène est présent dans le génome
				matrixByVfFamilies[familyName][genome] += matrixAllGenes[geneName][genome]  # incrémente la matrice de sa valeur

	return matrixByVfFamilies



def getResfinderMatrixByGenesTypes(matrixAllGenes, dicoGenes) :
	genomes = matrixAllGenes.index
	genesList = matrixAllGenes.columns
	antibioticFamilies = []

	for geneName in dicoGenes :
		antibioticFamily = dicoGenes[geneName].antibioticFamily

		if antibioticFamily not in antibioticFamilies :
			antibioticFamilies.append(antibioticFamily)

	matrixByAntibioticFamilies = pd.DataFrame(0, index = genomes, columns = antibioticFamilies)

	for genome in genomes :
		for geneName in dicoGenes :
			familyName = dicoGenes[geneName].antibioticFamily

			if matrixAllGenes[geneName][genome] != 0 : # si le gène est présent dans le génome
				matrixByAntibioticFamilies[familyName][genome] += matrixAllGenes[geneName][genome]  # incrémente la matrice de sa valeur

	return matrixByAntibioticFamilies


def writeMatrix(matrix, matrixName, dirMatrix, index) :
# Fonction qui écrit la matrice dans un fichier
	matrix.to_csv(dirMatrix + matrixName, sep='\t', index = index) # écriture de la matrice dans un fichier tsv



def getVfdbCorrespondanceTable(matrixAllGenes, dicoGenes) :
# Fonction qui réalise une table de correspondance en attribuant à chaque gène son type de résistance et le nombre de fois qu'il ets retrouvé tous les génomes confondu

	genes = []
	vfNames = []
	vfFamilies = []
	speciesNames = []
	number = []

	corTable = pd.DataFrame(columns = ['genes', 'species', 'VF names', 'family names', 'number']) # Définition d'un dataframe vide avec 3 colonnes

	genesList = matrixAllGenes.columns

	for geneName in genesList : # pour chaque colonne de la matrice
			genes.append(geneName) # ajout du gène à la liste
			number.append(sum(matrixAllGenes[geneName])) # exemplaire du gène dans tous les génomes

			vfName = dicoGenes[geneName].vfName
			vfNames.append(vfName)

			species = dicoGenes[geneName].species
			speciesNames.append(species)

			vfFamily = dicoGenes[geneName].vfFamily
			vfFamilies.append(vfFamily)

	corTable['genes'] = genes # remplissage de la colonne gene
	corTable['species'] = speciesNames
	corTable['VF names'] = vfNames # remplissage de la colonne antibiotique
	corTable['family names'] = vfFamilies
	corTable['number'] = number # remplissage de la colonne effectif

	print(corTable)
	print(corTable['number'].sum())

	return corTable

def getResfinderCorrespondanceTable(matrixAllGenes, dicoGenes) :
# Fonction qui réalise une table de correspondance en attribuant à chaque gène son type de résistance et le nombre de fois qu'il ets retrouvé tous les génomes confondu

	genes = []
	antibioticFamilies = []
	number = []

	corTable = pd.DataFrame(columns = ['genes', 'antibiotic resistance', 'number']) # Définition d'un dataframe vide avec 3 colonnes
	genesList = matrixAllGenes.columns

	for gene in genesList : # pour chaque colonne de la matrice

		if not gene in genes : # si le gène n'est pas dans la liste
			genes.append(gene) # ajout du gène à la liste
			number.append(sum(matrixAllGenes[gene])) # exemplaire du gène dans tous les génomes

			antibioticFamily = dicoGenes[gene].antibioticFamily
			antibioticFamilies.append(antibioticFamily)

	corTable['genes'] = genes # remplissage de la colonne gene
	corTable['antibiotic resistance'] = antibioticFamilies # remplissage de la colonne antibiotique
	corTable['number'] = number # remplissage de la colonne effectif

	print(corTable)
	print(corTable['number'].sum())

	return corTable


def getCorrespondanceTable(matrixAllGenes, dicoGenes) :
# Fonction qui réalise une table de correspondance en attribuant à chaque gène son type de résistance et le nombre de fois qu'il ets retrouvé tous les génomes confondu

	genes = []
	number = []

	corTable = pd.DataFrame(columns = ['genes', 'number']) # Définition d'un dataframe vide avec 3 colonnes
	
	genesList = matrixAllGenes.columns

	for gene in genesList : # pour chaque colonne de la matrice

		if not gene in genes : # si le gène n'est pas dans la liste
			genes.append(gene) # ajout du gène à la liste
			number.append(sum(matrixAllGenes[gene])) # exemplaire du gène dans tous les génomes

	corTable['genes'] = genes # remplissage de la colonne gene
	corTable['number'] = number # remplissage de la colonne effectif

	print(corTable)
	print(corTable['number'].sum())

	return corTable


def removeGenesPresentInAllGenomes(matrixAllGenes) :
# Fonction qui supprime dans la matrice les gènes présent dans tous les génomes et en fait le liste
	genesList = matrixAllGenes.columns # nom des gènes

	removedGenes = [] # liste des gènes supprimés

	for gene in genesList : # pour chaque gène
		print(gene)
		print(all(MATRIX[gene] == 1)) 
		if all(matrixAllGenes[gene] == 1) : # Si toutes les valuers de la colonne du gène sout égales à 1
			matrixAllGenes.pop(gene) # suppression du gène
			removedGenes.append(gene) # ajout du gène à la liste des gènes supprimés

	print(removedGenes)

	return matrixAllGenes


def heatmap(matrix, heatmapName, dirHeatmap) :
# fonction qui réalise une heatmap

	# sns.clustermap(df, metric="correlation", method="single", cmap="Blues", standard_scale=1, row_colors=row_colors)
	
	try : # si pas de RecursionError

		heatmap = sns.clustermap(matrix, metric="euclidean", method="average", figsize=(18, 14), linewidths=.003) # réalisation de la heatmap
		heatmap.savefig(dirHeatmap + heatmapName) # sauvegarde de la heatmap dans un png

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
	RESDIR = Arguments.resdir

	if WORKDIR[-1] != '/' :
		WORKDIR += '/'

	if RESDIR[-1] != '/' :
		RESDIR += '/'


	if not os.path.exists(WORKDIR + RESDIR) :
		os.system('mkdir ' + WORKDIR + RESDIR)
		os.system('mkdir ' + WORKDIR + RESDIR + 'matrix/')
		os.system('mkdir ' + WORKDIR + RESDIR + 'heatmap/')

	else : 
		os.system('mkdir ' + WORKDIR + RESDIR + 'matrix/')
		os.system('mkdir ' + WORKDIR + RESDIR + 'heatmap/')


	DIR_MATRIX = WORKDIR + RESDIR + "matrix/" # chemin de destination des matrices
	DIR_HEATMAP = WORKDIR + RESDIR + "heatmap/" #chemin de derstinatin des heatmap


	dicoGenomes = {}
	dicoGenes = {}

	if Arguments.defaultDatabase is not None :
		if Arguments.defaultDatabase == 'resfinder' : 
			getGenomesObjects(Arguments.abricateList, dicoGenomes, Arguments.defaultDatabase, dicoGenes, None)
			matrixAllGenes = getMatrixAllGenes(dicoGenomes)
			matrixByGenesTypes = getResfinderMatrixByGenesTypes(matrixAllGenes, dicoGenes)
			corTable = getResfinderCorrespondanceTable(matrixAllGenes, dicoGenes)

			writeMatrix(matrixByGenesTypes, "matrix_by_gene_type.tsv", DIR_MATRIX, True) # ecrit la matrice des antibio

			heatmap(matrixByGenesTypes, "heatmap_by_gene_type.png", DIR_HEATMAP) # réalise la heatmap des antibio

		elif Arguments.defaultDatabase == 'vfdb' : 
			dicoVfFamilies = getVfFamiliesDico('VFs.csv')
			print(dicoVfFamilies)
			getGenomesObjects(Arguments.abricateList, dicoGenomes, Arguments.defaultDatabase, dicoGenes, dicoVfFamilies)
			matrixAllGenes = getMatrixAllGenes(dicoGenomes)
			matrixByGenesTypes = getVfdbMatrixByGenesTypes(matrixAllGenes, dicoGenes)
			corTable = getVfdbCorrespondanceTable(matrixAllGenes, dicoGenes)

			writeMatrix(matrixByGenesTypes, "matrix_by_gene_type.tsv", DIR_MATRIX, True) # ecrit la matrice des antibio

			heatmap(matrixByGenesTypes, "heatmap_by_gene_type.png", DIR_HEATMAP) # réalise la heatmap des antibio

		else : 
			getGenomesObjects(Arguments.abricateList, dicoGenomes, Arguments.defaultDatabase, dicoGenes, None)
			matrixAllGenes = getMatrixAllGenes(dicoGenomes)
			corTable = getCorrespondanceTable(matrixAllGenes, dicoGenes)

		if Arguments.remove : 
			matrixAllGenes = removeGenesPresentInAllGenomes(matrixAllGenes)


	elif Arguments.privateDatabase : 
		getGenomesObjects(Arguments.abricateList, dicoGenomes, Arguments.defaultDatabase, dicoGenes, None)
		matrixAllGenes = getMatrixAllGenes(dicoGenomes)
		corTable = getCorrespondanceTable(matrixAllGenes, dicoGenes)

	writeMatrix(matrixAllGenes, "matrixAllGenes.tsv", DIR_MATRIX, True) # ecrit la matrice des gènes
	writeMatrix(corTable, "correspondance_table.tsv", DIR_MATRIX, False) # ecrit la table de correspondance

	heatmap(matrixAllGenes, "heatmap_all_genes.png", DIR_HEATMAP) # réalise la heatmap des gène

	#os.system("rm " + Arguments.abricateList)
		

	t2 = time.time()

	diff = round(t2 - t1,3)
	print ("Temps : " + str(diff) + " secondes")

# lancer la fonction main()  au lancement du script
if __name__ == "__main__":
	main()	  
