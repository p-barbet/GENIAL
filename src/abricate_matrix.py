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
							'ecoli_vf', 'plasmidfinder', 'vfdb', 'ncbi', 'vir_clost', 'enterotox_staph', 'phages'], help='default \
								database to use (resfinder, card, argannot, ecoh, ecoli_vf, plasmidfinder, vfdb, ncbi. Incompatible with --privatedb)')

	db_type.add_argument('--privatedb', action="store_true", dest='privateDatabase', \
						help='private database name. Implies -dbp, -dbf. Incompatible with --defaultdb')

	# parser.add_argument('--research', dest='research', choices=['antibiotic_resistance','virulence_factor'], required=True,
	# 				help='Type of researched genes antibiotic antibiotic resistance or virulence factor ')

	parser.add_argument('--R', action="store_true", dest='remove',
					default=False, help='remove genes present in all genomes from the matrix (default:False)')

	return parser


# objet génome (atributs : ID, fichier abricate, gènes, matrice abricate, matrice des gènes)
class genome(object) :

	def __init__(self) :
		self.ID = ""
		self.abricateFile = ""
		self.genes = [] 
		self.abricateMatrix = None # informations du fichier abricate
		self.genesMatrix = None


	def setID(self, ID) :
		self.ID = ID

	def setAbricateFile(self, abricateFile) :
		self.abricateFile = abricateFile

	def setAbricateMatrix(self) :
		self.abricateMatrix = pd.read_csv(self.abricateFile, sep='\t', index_col=0, dtype = str) # dataframe contenant les informations du fichie

	def setGenes(self) :

		genesList = self.abricateMatrix['GENE']

		if len(genesList) != 0 : # si au moins 1 gène de résistance trouvé
		
			for geneName in genesList :
				self.genes.append(geneName) # ajout de chaque gène à la liste gène

	def setGenesMatrix(self) :
		if len(self.genes) != 0 : # si au moins 1 gène de résistance trouvé
		
			for geneName in self.genes :
				self.genesMatrix = pd.DataFrame(1, index = [self.ID], columns = self.genes) # création d'un dataframe remplis de 1 avec les gènes du génome en colonnes et le nom du génome en index
				self.genesMatrix = self.genesMatrix.groupby(self.genesMatrix.columns, axis=1).sum() # Somme des colonnes avec le même nom de gene (prise en compte de gène trouvés plusieurs fois dans le génome)

		else :

			self.genesMatrix = pd.DataFrame(index = [self.ID]) # sinon création d'un dataframe vide avec l'ID du génome en index


	def getGenesObjects(self, dicoGenes, database, dicoVfFamilies) :
		for geneName in self.genes :
			if geneName not in dicoGenes : # si le nom du gène n'est pas dans le dictionnaire des gènes création de l'objet gène correspondant
				dicoGenes[geneName] = gene()

				

				dicoGenes[geneName].setName(geneName)
				dicoGenes[geneName].setDatabase(database)

				if database == 'vfdb' : 
					productInfos = self.abricateMatrix.set_index('GENE')['PRODUCT'][geneName] # récupération des informations dans la colonne produit du gène
					dicoGenes[geneName].setVfName(productInfos) # nom du facteur de virulence
					dicoGenes[geneName].setSpecies(productInfos) # espèce du facteur de virulence
					dicoGenes[geneName].setVfFamily(dicoVfFamilies) # famille du facteur de virulence

				elif database == 'resfinder' :
					dicoGenes[geneName].setAntibioticFamily(productInfos) # famille d'antibiotiques

	

# objet gène (atributs : nom, base de données, nom VF, nom famille VF, espèce, famille d'antibiotiques)
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
		if (self.species in dicoVfFamilies) and (self.vfName in dicoVfFamilies[self.species]) : # si l'espèce et facteur de vurlence sont dans le dictionnaire des facteurs de virulence
			self.vfFamily = dicoVfFamilies[self.species][self.vfName] # attribution de la famille du facteur de virulence

		else : 
			self.vfFamily = 'Unclassified'

	def setAntibioticFamily(self, productInfos) :
		self.antibioticFamily = productInfos



# Fonction qui crée tous les objets génomes et gènes et les stock dans des dictionnaires avec les IDs de ces derniers comme clés
def getGenomesObjects(inputFile, dicoGenomes, database, dicoGenes, dicoVfFamilies) :
	
	data = open(inputFile, 'r')
	lines = data.readlines() 
	data.close()

	for line in lines :

		line = line.rstrip() # suppression des retours chariots des lignes
		infos = line.split("\t") 

		ID = infos[1] # IDs des génomes
		abricateFile = infos[0] # chemins des résultats abricate

		dicoGenomes[ID] = genome() 

		dicoGenomes[ID].setID(ID) 
		dicoGenomes[ID].setAbricateFile(abricateFile)
		dicoGenomes[ID].setAbricateMatrix()
		dicoGenomes[ID].setGenes()
		dicoGenomes[ID].setGenesMatrix()


		dicoGenomes[ID].getGenesObjects(dicoGenes, database, dicoVfFamilies)


# Fonction qui stock dans un dictionnaire les familles des facteurs de virulence associées à chaque espèces
def getVfFamiliesDico(vfFamiliesFile) :

	data = open(vfFamiliesFile, 'r')
	lines  = data.readlines()
	data.close()

	dicoVfFamilies = {} # dictionnaire des familles de facteur de virulence

	for line in lines :

		line = line.rstrip() # suppresion des retours chariot
		infos = line.split('\t')

		vfName = infos[0] # nom du facteur de virulence
		species = infos[1] # nom de l'espèce

		try :
			vfFamily = infos[2].split('; ')[0] # famille du facteur de virulence

		except IndexError:
			vfFamily = 'Unclassified'

		if species in dicoVfFamilies : # si l'espèce est dans le dictionnaire
			dicoVfFamilies[species][vfName] = vfFamily # association de la famille au facteur de virulence

		else :
			dicoVfFamilies[species] = {} # sinon création de l'espèce
			dicoVfFamilies[species][vfName] = vfFamily # association de la famille au facteur de virulence

	return dicoVfFamilies



# Fonction qu réalise la matrice de présence abscence de tous les gènes pour chaque génome
def getMatrixAllGenes(dicoGenomes, dicoGenes, database) :

	dataframes = [] # liste des dataframes des génomes

	# construction de la liste des dataframes de tous les génomes
	for genome in dicoGenomes :
		dataframes.append(dicoGenomes[genome].genesMatrix) 

	matrixAllGenes = pd.concat(dataframes, sort=True) # concaténation des dataframes de la liste

	matrixAllGenes.fillna(0, inplace=True) # remplacement des NaN par des 0

	if database == "resfinder" or database == "vfdb" :
		matrixAllGenes = sortGenesByTypes(matrixAllGenes, dicoGenes, database) # tri des gènes par type si la base de donnée est resfinder ou vfdb

	return matrixAllGenes


# Fonction qui tri les gènes d'une matrice par type
def sortGenesByTypes(matrixAllGenes, dicoGenes, database) :

	genesList = matrixAllGenes.columns # liste des gènes
	newGenesList = [] # nouvelle liste des gènes réordonés
	

	# base de données resfinder
	if database == "resfinder" :
		antibioticFamilies = [] # liste des familles d'antibiotiques

		for geneName in genesList : 
			antibioticFamily = dicoGenes[geneName].antibioticFamily + "|" # famille d'antibiotique
			antibioticFamilies.append(antibioticFamily) # construction de la liste des familles d'antibiotiques

		genesAndFamilies = antibioticFamilies + genesList # association de chaque gènes à sa famille d'antibiotique


	# base de données vfdb
	elif database == "vfdb" : 
		vfFamilies = [] # liste des familles des facteurs de virulence

		for geneName in genesList : 
			vfFamily = dicoGenes[geneName].vfFamily + "|" + dicoGenes[geneName].vfName + "|" # famille et nom du facteur de virulence
			vfFamilies.append(vfFamily) # construction de la liste des familles des facteurs de virulence

		genesAndFamilies = vfFamilies + genesList # association de chaque gène à sa famille de facteur de virulence


	matrixAllGenes.columns = genesAndFamilies # renomage des gènes par ces derniers associés à leur famille

	matrixAllGenes = matrixAllGenes.sort_index(axis = 1) # tri des gènes par ordre alphabétique

	genesAndFamilies = matrixAllGenes.columns # gènes de la matrice

	for genesAndFamily in genesAndFamilies :
		geneName = genesAndFamily.split("|")[-1] # gènes sans sa famille
		newGenesList.append(geneName) # construction de la liste des gènes réordonnés

	matrixAllGenes.columns = newGenesList # renomage des gènes réordonés

	return matrixAllGenes



# Fonction qui construit la matrice par type de gènes pour la base de données vfdb
def getVfdbMatrixByGenesTypes(matrixAllGenes, dicoGenes) :

	genomes = matrixAllGenes.index # génomes
	genesList = matrixAllGenes.columns # gènes

	vfFamilies = [] # liste des familles des facteurs de virulence

	for geneName in genesList : # pour chaqua gène
		vfFamily = dicoGenes[geneName].vfFamily # famille du gène

		if vfFamily not in vfFamilies :
			vfFamilies.append(vfFamily)  # construction de la liste des familles des facteurs de virulence

	matrixByVfFamilies = pd.DataFrame(0, index = genomes, columns = vfFamilies) # construction d'une matrice remplie de 0 avec les génomes en index et les familles en colonne

	for genome in genomes : # pour chaque génome
		for geneName in genesList :
			familyName = dicoGenes[geneName].vfFamily  # famille du gène

			if matrixAllGenes[geneName][genome] != 0 : # si le gène est présent dans le génome
				matrixByVfFamilies[familyName][genome] += matrixAllGenes[geneName][genome]  # incrémentation de la matrice du nombre de fois qu'on le retrouve dans le génome

	return matrixByVfFamilies


# Fonction qui construit la matrice par type de gènes pour la base de données resfinder
def getResfinderMatrixByGenesTypes(matrixAllGenes, dicoGenes) :

	genomes = matrixAllGenes.index # genomes
	genesList = matrixAllGenes.columns # gènes

	antibioticFamilies = [] # liste des familles d'antibiotiques

	for geneName in dicoGenes : # pour chaque gène
		antibioticFamily = dicoGenes[geneName].antibioticFamily # famille d'antibiotique

		if antibioticFamily not in antibioticFamilies : 
			antibioticFamilies.append(antibioticFamily)# construction de la liste des familles d'antibiotiques

	matrixByAntibioticFamilies = pd.DataFrame(0, index = genomes, columns = antibioticFamilies) # construction d'une matrice remplie de 0 avec les génomes en index et les familles en colonne

	for genome in genomes :
		for geneName in dicoGenes :
			familyName = dicoGenes[geneName].antibioticFamily # femille d'antibiotique

			if matrixAllGenes[geneName][genome] != 0 : # si le gène est présent dans le génome
				matrixByAntibioticFamilies[familyName][genome] += matrixAllGenes[geneName][genome] # incrémentation de la matrice du nombre de fois qu'on le retrouve dans le génome

	return matrixByAntibioticFamilies



# Fonction qui écrit une latrice dans un fichier tabulé (tsv)
def writeMatrix(matrix, matrixName, dirMatrix, index) :
	matrix.to_csv(dirMatrix + matrixName, sep='\t', index = index)


# Fonction qui construit une table de correspondance en attribuant à chaque gène son son espèces, son VF, sa famile et son effectif pour la base de données vfdb
def getVfdbCorrespondanceTable(matrixAllGenes, dicoGenes) :

	genes = [] # gènes
	vfNames = [] # noms des VF
	vfFamilies = [] # noms des familles des VF
	speciesNames = [] # nom des espèces
	number = [] # effectifs des gènes

	corTable = pd.DataFrame(columns = ['genes', 'species', 'VF names', 'family names', 'number']) # dataframe vide avec 5 colonnes

	genesList = matrixAllGenes.columns # gènes

	for geneName in genesList : # pour chaque gène
			genes.append(geneName) # ajout du gène à la liste des gènes

			nb = sum(matrixAllGenes[geneName]) # effectif du gène
			number.append(nb) # ajout à la liste

			vfName = dicoGenes[geneName].vfName # nom du VF
			vfNames.append(vfName) # ajout à la liste

			species = dicoGenes[geneName].species # nom de l'espèce
			speciesNames.append(species) # ajout à la liste

			vfFamily = dicoGenes[geneName].vfFamily # nom d ela famille du VF
			vfFamilies.append(vfFamily) # ajout à la liste

	corTable['genes'] = genes # remplissage de la colonne gene
	corTable['species'] = speciesNames # remplissage de la colonne espèce
	corTable['VF names'] = vfNames # remplissage de la colonne antibiotique
	corTable['family names'] = vfFamilies # remplissage de la colonne famille du VF
	corTable['number'] = number # remplissage de la colonne effectif

	print(corTable)
	print(corTable['number'].sum())

	return corTable


# Fonction qui construit une table de correspondance en attribuant à chaque gène son type et son effectif pour la base de données resfinder
def getResfinderCorrespondanceTable(matrixAllGenes, dicoGenes) :

	genes = [] # gènes
	antibioticFamilies = [] # familles d'natibiotiques
	number = [] # effecyifs des gènes

	corTable = pd.DataFrame(columns = ['genes', 'antibiotic resistance', 'number']) # dataframe vide avec 3 colonnes
	genesList = matrixAllGenes.columns

	for geneName in genesList : # pour chaque gène

		if not geneName in genes : 
			genes.append(geneName) # ajout du gène à la liste

			nb = sum(matrixAllGenes[geneName]) # effectif du gène
			number.append(nb) # ajout à la liste

			antibioticFamily = dicoGenes[geneName].antibioticFamily # famille du gène
			antibioticFamilies.append(antibioticFamily) # ajout à la liste

	corTable['genes'] = genes # remplissage de la colonne gene
	corTable['antibiotic resistance'] = antibioticFamilies # remplissage de la colonne antibiotique
	corTable['number'] = number # remplissage de la colonne effectif

	print(corTable)
	print(corTable['number'].sum())

	return corTable



# Fonction qui construit une table de correspondance en attribuant à chaque gène son effectif pour la base de données resfinder
def getCorrespondanceTable(matrixAllGenes, dicoGenes) :

	genes = [] # genes
	number = [] # effectifs des gènes

	corTable = pd.DataFrame(columns = ['genes', 'number']) # dataframe vide avec 2 colonnes
	
	genesList = matrixAllGenes.columns # liste des gènes

	for geneName in genesList : # pour chaque gène

		if not geneName in genes :
			genes.append(geneName) # ajout du gène à la liste

			nb = sum(matrixAllGenes[geneName]) # effectif du gène
			number.append(nb) # ajout à la liste


	corTable['genes'] = genes # remplissage de la colonne gene
	corTable['number'] = number # remplissage de la colonne effectif

	print(corTable)
	print(corTable['number'].sum())

	return corTable



# Fonction qui supprime les gènes présents dans tous les génomes 
def removeGenesPresentInAllGenomes(matrixAllGenes) :

	genesList = matrixAllGenes.columns # gènes

	removedGenes = [] # liste des gènes supprimés

	for geneName in genesList : # pour chaque gène 
		if all(matrixAllGenes[geneName] == 1) : # Si toutes les valeurs de la colonne du gène sout égales à 1
			matrixAllGenes.pop(geneName) # suppression du gène
			removedGenes.append(geneName) # ajout du gène à la liste des gènes supprimés

	return matrixAllGenes



# Fonction qui construit une heatmap à partir d'une heatmap
def heatmap(matrix, heatmapName, dirHeatmap, nbGenomes) :
	
	try : 
		heatmap = sns.clustermap(matrix, metric="euclidean", method="average", figsize=(round(0.23*nbGenomes + 0.66), round(0.23*nbGenomes + 0.66)), linewidths=.003, col_cluster = False, yticklabels=True, cmap = "Greys")  # construction de la heatmap
		heatmap.savefig(dirHeatmap + heatmapName) # sauvegarde de la heatmap dans un png

	except RecursionError: 
		print("Erreur lors de la construction de la heatmap, la matrice est trop grande") # affichage d'un message d'erreur si la matrice est trop grande



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
			getGenomesObjects(Arguments.abricateList, dicoGenomes, Arguments.defaultDatabase, dicoGenes, None) # construction des objets génomes
			matrixAllGenes = getMatrixAllGenes(dicoGenomes, dicoGenes, Arguments.defaultDatabase) # construction de la matrice avec tous les gènes
			matrixByGenesTypes = getResfinderMatrixByGenesTypes(matrixAllGenes, dicoGenes) # construction de la matrica par type de gène
			corTable = getResfinderCorrespondanceTable(matrixAllGenes, dicoGenes) # construction de la table de correspondance

			writeMatrix(matrixByGenesTypes, "matrix_by_gene_type.tsv", DIR_MATRIX, True) # ecriture de la matrice par type de gènes

			heatmap(matrixByGenesTypes, "heatmap_by_gene_type.png", DIR_HEATMAP, len(dicoGenomes)) # construction de la heatmap par type de gènes

		elif Arguments.defaultDatabase == 'vfdb' : 
			dicoVfFamilies = getVfFamiliesDico('VFs.tsv') # construction du dictionnaire des familles de VF
			getGenomesObjects(Arguments.abricateList, dicoGenomes, Arguments.defaultDatabase, dicoGenes, dicoVfFamilies) # construction des objets génomes
			matrixAllGenes = getMatrixAllGenes(dicoGenomes, dicoGenes, Arguments.defaultDatabase) # construction de la matrice avec tous les gènes
			matrixByGenesTypes = getVfdbMatrixByGenesTypes(matrixAllGenes, dicoGenes) # construction de la matrica par type de gène
			corTable = getVfdbCorrespondanceTable(matrixAllGenes, dicoGenes) # construction de la table de correspondance

			writeMatrix(matrixByGenesTypes, "matrix_by_gene_type.tsv", DIR_MATRIX, True) # ecriture de la matrice par type de gènes

			heatmap(matrixByGenesTypes, "heatmap_by_gene_type.png", DIR_HEATMAP, len(dicoGenomes)) # construction de la heatmap par type de gènes

		else : 
			getGenomesObjects(Arguments.abricateList, dicoGenomes, Arguments.defaultDatabase, dicoGenes, None) # construction des objets génomes
			matrixAllGenes = getMatrixAllGenes(dicoGenomes, dicoGenes, Arguments.defaultDatabase) # construction de la matrice avec tous les gènes
			corTable = getCorrespondanceTable(matrixAllGenes, dicoGenes) # construction de la table de correspondance

		if Arguments.remove : 
			matrixAllGenes = removeGenesPresentInAllGenomes(matrixAllGenes) # retrait des gènes présent dans tous les génomes


	elif Arguments.privateDatabase : 
		getGenomesObjects(Arguments.abricateList, dicoGenomes, Arguments.defaultDatabase, dicoGenes, None) # construction des objets génomes
		matrixAllGenes = getMatrixAllGenes(dicoGenomes, dicoGenes, Arguments.defaultDatabase) # construction de la matrice avec tous les gènes
		corTable = getCorrespondanceTable(matrixAllGenes, dicoGenes) # construction de la table de correspondance

	writeMatrix(matrixAllGenes, "matrixAllGenes.tsv", DIR_MATRIX, True) # ecriture de la matrice pour tous les gènes
	writeMatrix(corTable, "correspondance_table.tsv", DIR_MATRIX, False) # ecriture la table de correspondance

	heatmap(matrixAllGenes, "heatmap_all_genes.png", DIR_HEATMAP, len(dicoGenomes)) # construction de la heatmap pour tous les gènes

	print(matrixAllGenes)

	os.system("rm " + Arguments.abricateList)
		

	t2 = time.time()

	diff = round(t2 - t1,3)
	print ("Temps : " + str(diff) + " secondes")

# lancer la fonction main()  au lancement du script
if __name__ == "__main__":
	main()	  

