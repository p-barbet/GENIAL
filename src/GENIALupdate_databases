#!/usr/bin/python3
# -*- coding: iso-8859-1 -*-

import os, sys, time
import pandas as pd
import argparse
from Bio import SeqIO


def get_parser() :
	# Fonction permettant de pourvoir demander des arguments

	parser = argparse.ArgumentParser(description= \
		"Update vfdb or resfinder databases", \
			formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	parser.add_argument("-db", action="store", dest="database", required = True,\
						type=str, choices=["resfinder", "vfdb"], help="database to be updated")

	parser.add_argument("-w", action="store", dest="workdir", 
					type=str, default=".", help="working directory")

	parser.add_argument("-dbp", action="store", dest="dbPath", 
					type=str, required=True, help="abricate databases repertory")

	return parser

# Objet séquence (atributs : header, sequence, newheader)
class Sequence(object) :

	def __init__(self) :
		self.header = ""
		self.sequence = []
		self.newHeader = ""

	def setHeader(self, header) :
		self.header = header

	def setSequence(self, sequence) :
		self.sequence = sequence

	def setNewHeader(self, newHeader) :
		self.newHeader = newHeader


# Fonction qui télécharge un fichier
def downloadFile(url, resdir, database) :

	compressedFileName = resdir + url.split("/")[-1]

	os.system("wget -O " + compressedFileName + " " + url)

	if database == "vfdb" :
		os.system("gunzip " + compressedFileName)
		
	else :
		os.system("unzip " + compressedFileName)
		os.system("rm " + compressedFileName)


# Fonction qui met à jour le fichier des facteurs de virulence
def vfFileUpdate(resdir) :

	vfFileUrl = "www.mgc.ac.cn/VFs/Down/VFs.xls.gz"
	downloadFile(vfFileUrl, resdir, "vfdb")
	FileName = resdir + "VFs.xls"
	excelFile = pd.read_excel(FileName, header = 1)
	columns = excelFile.columns

	for column in columns : 

		if column not in ["VF_Name", "Bacteria", "Keyword"] :
			excelFile.pop(column)

	csvFile = resdir + "VFs.tsv"
	excelFile.to_csv(csvFile, sep="\t", index=False)
	os.system("rm " + FileName)

	

# Fonction qui met à jour vfdb
def getVfdbSequences(resdir, dicoSequences) :

	vfdbFastaUrl = "www.mgc.ac.cn/VFs/Down/VFDB_setA_nt.fas.gz"
	downloadFile(vfdbFastaUrl, resdir, "vfdb")

	vfdbFastaName = resdir + "VFDB_setA_nt.fas"

	for Seq in SeqIO.parse(vfdbFastaName, "fasta") :

		header = Seq.description
		geneName = header.split(" ")[1].replace("(", "").replace(")", "")

		try :
			accession = header.split(" ")[0].split("|")[1].replace(")", "")

		except IndexError:
			accession = ""

		info = " ".join(header.split(" ")[1:])

		newHeader = ">vfdb" + "~~~" + geneName + "~~~" + accession + " " + info

		sequence = str(Seq.seq)
		
		if header not in dicoSequences :
			dicoSequences[header] = Sequence()

			dicoSequences[header].setHeader(header)
			dicoSequences[header].setNewHeader(newHeader)
			dicoSequences[header].setSequence(sequence)

		else :
			if dicoSequences[header].sequence == sequence :
				break

			else :
				print("\n" + newHeader)

	os.system("rm " + vfdbFastaName)



def getResfinderSequences(resdir, dicoSequences) :

	resfinderFastaUrl = "https://bitbucket.org/genomicepidemiology/resfinder_db/get/d3d7a6ceaa49.zip"
	downloadFile(resfinderFastaUrl, resdir, "resfinder")
	os.system("mv genomicepidemiology-resfinder_db-d3d7a6ceaa49/ " + resdir)

	resfinderDirectoryName = resdir + "genomicepidemiology-resfinder_db-d3d7a6ceaa49/"

	antibioticsList = ["aminoglycoside", "beta-lactam", "colistin", "fosfomycin", "fusidicacid", "glycopeptide", "macrolide", "nitroimidazole", "oxazolidinone", "phenicol", "quinolone", "rifampicin", "sulphonamide", "tetracycline", "trimethoprim"]

	for antibiotic in antibioticsList :

		for Seq in SeqIO.parse(resfinderDirectoryName + antibiotic + ".fsa", "fasta") :

			header = Seq.description
			geneName = "_".join(header.split("_")[0:-1])

			try :
				accession = header.split("_")[-1]

			except IndexError:
				accession = ""

			info = antibiotic

			dicoCle = ">resfinder" + "~~~" + geneName + "~~~" + accession

			newHeader = dicoCle + " " + info

			sequence = str(Seq.seq)

			if dicoCle not in dicoSequences :
				dicoSequences[dicoCle] = Sequence()

				dicoSequences[dicoCle].setHeader(header)
				dicoSequences[dicoCle].setNewHeader(newHeader)
				dicoSequences[dicoCle].setSequence(sequence)

			else :
				if dicoSequences[dicoCle].sequence == sequence and dicoSequences[dicoCle].newHeader.split(" ")[-1] != info :
					dicoSequences[dicoCle].newHeader = dicoSequences[dicoCle].newHeader + "," + info

				elif dicoSequences[dicoCle].sequence == sequence :
					break

				else :
					print("\n" + newHeader)

					
	os.system("rm -R " + resdir + "genomicepidemiology-resfinder_db-d3d7a6ceaa49/")


def writeFinalFasta(dicoSequences, resdir) :
	newVfdbFasta = open(resdir + "sequences", "w")

	for sequence in dicoSequences :
		newVfdbFasta.write(dicoSequences[sequence].newHeader + "\n")

		sequenceLength = len(dicoSequences[sequence].sequence)

		sequencePosition = 0

		while sequencePosition < sequenceLength :

			newVfdbFasta.write(dicoSequences[sequence].sequence[sequencePosition:sequencePosition+60] + "\n")
			sequencePosition += 60

	newVfdbFasta.close()


def setupdb(dbName, abricateDbsRepertory, dbMultifasta) :

	multifastaName = dbMultifasta.split("/")[-1] # nom du fichier multifasta de la base a créer

	if abricateDbsRepertory[-1] != "/" :
		abricateDbsRepertory += "/"

	if not os.path.exists(abricateDbsRepertory + dbName) : # création du répertoire de la base dans abricate si il n"existe pas
		os.system("mkdir " + abricateDbsRepertory + dbName) 

	os.system("cp " + dbMultifasta + " " + abricateDbsRepertory + dbName) # copie du multifasta dans le répartoire abricate

	os.system("abricate --setupdb") # idexation de la base dans abricate


def main():
	
	##################### gets arguments #####################

	parser=get_parser()
	
	#print parser.help if no arguments
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	
	# mettre tout les arguments dans la variable Argument
	Arguments=parser.parse_args()

	WORKDIR = Arguments.workdir

	if WORKDIR[-1] != "/" :
		WORKDIR += "/"

	SEQUENCES = {}

	DATABASE = Arguments.database
	RESDIR = WORKDIR + DATABASE + "_update/"

	if not os.path.exists(RESDIR) :
		os.system("mkdir " + RESDIR)

	dicoSequences = {}
	
	if DATABASE == "vfdb" :
		vfFileUpdate(RESDIR)
		getVfdbSequences(RESDIR, dicoSequences)
		writeFinalFasta(dicoSequences, RESDIR)
		setupdb(DATABASE, Arguments.dbPath, RESDIR + "sequences")

	if DATABASE == "resfinder" :
		getResfinderSequences(RESDIR, dicoSequences)
		writeFinalFasta(dicoSequences, RESDIR)
		setupdb(DATABASE, Arguments.dbPath, RESDIR + "sequences")


# lancer la fonction main()  au lancement du script
if __name__ == "__main__":
	main()	  

