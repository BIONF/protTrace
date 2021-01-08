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
import argparse
import configure
import distanceCalculation
import preprocessing
import traceabilityCalculation
import mapToSpeciesTree
import time

def main(argv):

    def Argparse ():
        """ Parses the arguments into an object """

        parser = argparse.ArgumentParser(description='Calculates the evolutionary traceability of query proteins to species specified in the species mapping file.')

        mandatory_arguments = parser.add_argument_group(title='Required arguments')
        input_group = mandatory_arguments.add_mutually_exclusive_group(required=True)
        input_group.add_argument('-i','--id',type=str,help='Path to a text file containing protein IDs')
        input_group.add_argument('-f','--fasta',type=str,help='Path to a text file containing proteins in FASTA format')  

        mandatory_arguments.add_argument('-c','--config',type=str,required=True,help='Path to the configuration file which can be created with the bin/create_config.pl script')

        modes = parser.add_argument_group(title='Additional modes',description='These options are used for maintenance.')      
        modes.add_argument('-d','--distance',action='store_true',help='Only updates the ML distances between all species in the species mapping file')

        return parser.parse_args()

    arguments = Argparse()

    id_list, fasta_list = '', ''
    only_update_all_distances = False

    # One of these arguments is required already
    if arguments.id:
        id_list = arguments.id
    if arguments.fasta:
        fasta_list = arguments.fasta

    # The config argument is required
    config_file = arguments.config

    if arguments.distance:
        only_update_all_distances = True

    config_file = os.path.abspath(config_file)

    # Calling the class in configure.py module and setting the tool parameters
    proteinParams = configure.setParams(config_file)

    sys.exit()

    # This is a special setting, where no protein ID is needed and its only
    # purpose is to update the distances between all species in the species list
    if only_update_all_distances:
        # This argument hopefully breaks anything but the intended routine
        proteinParams.species = 'ALL'
        distanceCalculation.calculate_species_distances(proteinParams)
    else:    
        # Check available pairwise species distances and calculate missing ones
        distanceCalculation.calculate_species_distances(proteinParams)
    
        if id_list != '':
            with open(id_list, 'r') as id_file:
                for line in id_file:
                    input_id = line.split()[0]
                    print('##### Running for OMA id: {0} #####'.format(input_id))
                    if proteinParams.preprocessing:
                        preprocessing.Preprocessing(input_id, 'None', config_file)
                    if proteinParams.traceability_calculation:
                        traceabilityCalculation.main(input_id, config_file)
                    if proteinParams.mapTraceabilitySpeciesTree:
                        mapToSpeciesTree.main(input_id, config_file)
        elif fasta_list != '':
            with open(fasta_list) as fa:
                for seqs in fa:
                    if '>' in seqs:
                        print('##### Running for fasta id: {0} #####'.format(seqs[1:-1]))
                        input_id = seqs.split()[0][1:]
                        query_seq = next(fa)
    
                    if proteinParams.preprocessing:
                        preprocessing.Preprocessing(input_id, query_seq, config_file)
                    if proteinParams.traceability_calculation:
                        traceabilityCalculation.main(input_id, config_file)
                    if proteinParams.mapTraceabilitySpeciesTree:
                        mapToSpeciesTree.main(input_id, config_file)

if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        print('ERROR:\tNo arguments entered for the traceability run:\nUSAGE:\tprotTrace.py -i <omaIdsFile> | -f <fastaSeqsFile> -c <configFile> [-help]')
        sys.exit(2)
    else:
        start_time = time.time()
        print('##### Start time: %s #####' %time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.localtime()))
        main(sys.argv[1:])
        end_time = time.time()
        print('##### End time: %s #####' %time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.localtime()))
        print('##### TOTAL TIME: %s hours#####' %((end_time - start_time) / 3600))
