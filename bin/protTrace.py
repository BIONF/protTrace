#!/usr/bin/env python

# Protein traceability prediction script (Main / Application)
# Author: Arpit Jain
# Date: 25 May, 2015
#
# Date: 01 September, 2015
# Modified the script to run more robustly with sequences file as input
#
# Date: 10 November, 2016
# HaMStR-OneSeq run bug fixed
# When performing HaMStR search, if we are left with just 1 sequence, automatically redirect to HaMStR-OneSeq search
# Sometimes input files had '\r' character before a new line ('\n'), this has been fixed while reading the input file
#

import os, sys
import getopt
import configure
import preprocessing
import traceabilityCalculation
import mapToSpeciesTree
import time

def main(argv):
	
	id_list, fasta_list, config_file = '', '', ''

	# Setting the get options method to read the input arguments
	try:
		opts, args = getopt.getopt(argv, "f:i:c:h", ["fasta=", "id=", "config=", "help"])
	except getopt.GetoptError:
		print 'Invalid arguments:\nUsage:\tprotTrace.py -i <omaIdsFile> | -f <fastaSeqsFile> -c <configFile> [-help]'
		sys.exit(2)

	for opt, arg in opts:
		if opt in ('-h','--help'):
			print "USAGE:\tprotTrace.py -i <omaIdsFile> | -f <fastaSeqsFile> -c <configFile> [-h]\n\t-i\t\tText file containing protein OMA ids (1 id per line)\n\t-f\t\tList of input protein sequences in fasta format\n\t-c\t\tConfiguration file for setting program's dependencies"
			sys.exit(2)
		elif opt in ('-i', '--id'):
			id_list = arg
		elif opt in ('-f','--fasta'):
			fasta_list = arg
		elif opt in ('-c','--config'):
			config_file = arg
		else:
			print 'Invalid arguments:\nUsage:\tprotTrace.py -i <omaIdsFile> | -f <fastaSeqsFile> -c <configFile> [-help]'
			sys.exit(2)
	
	config_file = os.path.abspath(config_file)
	
	# Calling the class in configure.py module and setting the tool parameters
	proteinParams = configure.setParams(config_file)
	
	if id_list != '':
		for ids in open(id_list):
			print '##### Running for OMA id: %s #####' %ids.split()[0]
			if proteinParams.preprocessing:	
				preprocessing.Preprocessing(ids.split()[0], 'None', config_file)
			if proteinParams.traceability_calculation:
				traceabilityCalculation.main(ids.split()[0], config_file)
			if proteinParams.mapTraceabilitySpeciesTree:
				mapToSpeciesTree.main(ids.split()[0], config_file)
	elif fasta_list != '':
		with open(fasta_list) as fa:
			for seqs in fa:
				if '>' in seqs:
					print '##### Running for fasta id: %s #####' %seqs[1:-1]
					inputId = seqs.split()[0][1:]
					querySeq = fa.next()

				if proteinParams.preprocessing:
					preprocessing.Preprocessing(inputId, querySeq, config_file)
				if proteinParams.traceability_calculation:
					traceabilityCalculation.main(inputId, config_file)
				if proteinParams.mapTraceabilitySpeciesTree:
					mapToSpeciesTree.main(inputId, config_file)
	

if __name__ == "__main__":
	if len(sys.argv[1:]) == 0:
		print 'ERROR:\tNo arguments entered for the traceability run:\nUSAGE:\tprotTrace.py -i <omaIdsFile> | -f <fastaSeqsFile> -c <configFile> [-help]'
		sys.exit(2)
	else:
		start_time = time.time()
		print '##### Start time: %s #####' %time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.localtime())
		main(sys.argv[1:])
		end_time = time.time()
		print '##### End time: %s #####' %time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.localtime())
		print '##### TOTAL TIME: %s hours#####' %((end_time - start_time) / 3600)
