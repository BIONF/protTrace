#!/usr/bin/env python
# -*- coding: utf-8 -*-

#######################################################################
# Copyright (C) 2020 Arpit Jain, Dominik Perisa,
# Prof. Dr. Ingo Ebersberger
#
#  This script is part of ProtTrace.
#
#  This script is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License <http://www.gnu.org/licenses/> for
#  more details
#
#  Contact: ebersberger@bio.uni-frankfurt.de
#
#######################################################################

# Protein traceability prediction script (Main / Application)
# Author: Arpit Jain
# Date: 25 May, 2015
#
# Date: 01 September, 2015
# Modified the script to run more robustly with sequences file as input
#
# Date: 10 November, 2016
# HaMStR-OneSeq run bug fixed
# When performing HaMStR search, if we are left with just 1 sequence,
# automatically redirect to HaMStR-OneSeq search. Sometimes input files had
# '\r' character before a new line ('\n'), this has been fixed while reading
# the input file
#

import os
import sys
import argparse
import configure
import distanceCalculation
import preprocessing
import traceabilityCalculation
import mapToSpeciesTree
from log import print_progress, time_report


def main(argv):

    def interpret_arguments():
        """ Parses the arguments into an object """

        parser = argparse.ArgumentParser(description=''
                                         'Calculates the evolutionary '
                                         'traceability of query proteins to '
                                         'species specified in the species '
                                         'mapping file.')

        mandatory_arguments = parser.add_argument_group(title=''
                                                        'Required arguments')
        input_group = mandatory_arguments.add_mutually_exclusive_group(
            required=True)
        input_group.add_argument('-i', '--id', type=str, help=''
                                 'Path to a text file containing protein IDs')
        input_group.add_argument('-f', '--fasta', type=str, help=''
                                 'Path to a text file containing proteins in '
                                 'FASTA format')

        mandatory_arguments.add_argument('-c', '--config', type=str,
                                         required=True, help='Path to the '
                                         'configuration file which can be '
                                         'created with the '
                                         'bin/create_config.pl script')

        modes = parser.add_argument_group(title='Additional modes',
                                          description='These options are '
                                          'used for maintenance.')
        modes.add_argument('-d', '--distance', action='store_true', help=''
                           'Only updates thie ML distances between all '
                           'species in the species mapping file')
        modes.add_argument('--ignore-distances', action='store_true', help=''
                           'Sets all pairwise species distances to 1, useful'
                           'for debugging the traceability computation '
                           'immediately.')

        return parser.parse_args()

    arguments = interpret_arguments()

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

    # This is a special setting, where no protein ID is needed and its only
    # purpose is to update the distances between all species in the species
    # list
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
                    print_progress('Running for OMA id {0}'.format(input_id))
                    if proteinParams.preprocessing:
                        preprocessing.Preprocessing(input_id, 'None',
                                                    config_file)
                    if proteinParams.traceability_calculation:
                        traceabilityCalculation.main(input_id, config_file)
                    if proteinParams.mapTraceabilitySpeciesTree:
                        mapToSpeciesTree.main(input_id, config_file)
        elif fasta_list != '':
            with open(fasta_list) as fa:
                for seqs in fa:
                    if '>' in seqs:
                        print_progress('Running for fasta id: {0}'
                                       .format(seqs[1:-1]))
                        input_id = seqs.split()[0][1:]
                        query_seq = next(fa)

                    if proteinParams.preprocessing:
                        preprocessing.Preprocessing(input_id, query_seq,
                                                    config_file)
                    if proteinParams.traceability_calculation:
                        traceabilityCalculation.main(input_id, config_file)
                    if proteinParams.mapTraceabilitySpeciesTree:
                        mapToSpeciesTree.main(input_id, config_file)


if __name__ == "__main__":
    if len(sys.argv[1:]) > 0:
        time_report.format_verbose('Start time')
        time_record = time_report()

        # Execute the software.
        main(sys.argv[1:])

        time_record.print_time('TOTAL TIME', 'hrs')
        time_report.format_verbose('End time')
