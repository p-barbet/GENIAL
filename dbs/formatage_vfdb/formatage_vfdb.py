#!/usr/bin/python3
# -*- coding: iso-8859-1 -*-

import os, sys, time


db_file = open('VFDB.fas', 'r')
db = db_file.readlines()
db_file.close()

new_db = open('sequences', 'w')

db_name = "vfdb2"

for line in db :
	line = line.rstrip()
	if line != '' and line[0] == '>' : 

		gene = line.split(' ')[1].replace('(', '').replace(')', '')

		try :
			accession = line.split(' ')[0].split('|')[1].replace(')', '')

		except IndexError:
			accession = ''

		info = ' '.join(line.split(' ')[1:])

		new_header = '>' + db_name + '~~~' + gene + '~~~' + accession + ' ' + info

		new_db.write(new_header + '\n')

	else : 
		new_db.write(line + '\n')

new_db.close()



 
