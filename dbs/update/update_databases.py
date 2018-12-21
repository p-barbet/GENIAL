#!/usr/bin/python3
# -*- coding: iso-8859-1 -*-

import os, sys, time
import pandas as pd
import argparse


def get_parser() :
	# Fonction permettant de pourvoir demander des arguments

	parser = argparse.ArgumentParser(description= \
		'Update vfdb database and virulence factors description file', \
			formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	parser.add_argument('-db', action="store", dest='database', required = True,\
						type=str, choices=['resfinder', 'vfdb'], help='database to be updated')


	parser.add_argument('-w', action="store", dest='workdir', 
					type=str, required=True, help='working directory')

	return parser


class sequence(object) :

	def __init__(self) :
		self.header = ""
		self.sequence = []
		self.geneName = ""
		self.accession = ""


	def setHeader(self, header) :
		self.header = header


	def setSequence(self, sequence) :
		self.sequence



def downloadFile(url, resdir) :

	compressedFileName = resdir + url.split('/')[-1]
	
	FileName = resdir + '.'.join(url.split('/')[-1].split('.')[0:2])

	os.system("wget -O " + compressedFileName + " " + url)

	os.system("gunzip -c " + compressedFileName + " > " + FileName)

	os.system("rm " + compressedFileName)

	return FileName



def vfFileUpdate(resdir) :

	vfFileUrl = "www.mgc.ac.cn/VFs/Down/VFs.xls.gz"

	vfFileName = downloadFile(vfFileUrl, resdir)
	
	excelFile = pd.read_excel(vfFileName, header = 1)

	columns = excelFile.columns

	for column in columns : 
		if column not in ["VF_Name", "Bacteria", "Keyword"] :
			excelFile.pop(column)

	csvFile = resdir + "VFs.tsv"

	excelFile.to_csv(csvFile, sep="\t", index=False)

	print(excelFile)


def vfdbDatabaseUpdate(resdir) :

	vfdbFastaUrl = "www.mgc.ac.cn/VFs/Down/VFDB_setA_nt.fas.gz"

	vfdbFastaName = downloadFile(vfdbFastaUrl, resdir)

	vfdbFastaFile = open(vfdbFastaName, 'r')
	vfdbFasta = vfdbFastaFile.readlines()
	vfdbFastaFile.close()

	newVfdbFasta = open(resdir + 'sequences', 'w')

	for line in vfdbFasta :
		line = line.rstrip()
		if line != '' and line[0] == '>' : 

			geneName = line.split(' ')[1].replace('(', '').replace(')', '')

			try :
				accession = line.split(' ')[0].split('|')[1].replace(')', '')

			except IndexError:
				accession = ''

			info = ' '.join(line.split(' ')[1:])

			newHeader = '> vfdb' + '~~~' + geneName + '~~~' + accession + ' ' + info

			newVfdbFasta.write(newHeader + '\n')

		else : 
			newVfdbFasta.write(line + '\n')
			
	newVfdbFasta.close()


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

	if WORKDIR[-1] != '/' :
		WORKDIR += '/'

	SEQUENCES = {}

	DATABASE = Argument.database
	RESDIR = WORKDIR + Database + "Update/"

	if not os.path.exists(RESDIR) :

		os.system('mkdir ' + RESDIR)
	
	if DATABASE == "vfdb" :

		vfFileUpdate(RESDIR)
		vfdbDatabaseUpdate(RESDIR, SEQUENCES)



# lancer la fonction main()  au lancement du script
if __name__ == "__main__":
	main()	  

