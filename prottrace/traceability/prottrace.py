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

""" The main entry point for ProtTrace which dispatches the query onto the
    different modules. """

from argparse import ArgumentParser
from pathlib import Path

from utils.configure import set_params
from utils.data_api import (
    prot_id,
    species_mapping,
    fasta,
    proteome
)
from utils.log import (
    print_progress,
    print_error,
    time_report
)

from species_distance.distance import calculate_species_distances

from traceability.evolutionary_model import main as evol_model
from traceability.self_hit_statistics import main as trace
from traceability.calculate_traceability import main as calc_trace
from traceability.colorize_tree import main as colorize_tree


def main():

    def interpret_arguments():
        """ Parses the arguments into an object """

        parser = ArgumentParser(description=''
                                'Calculates the evolutionary '
                                'traceability of query proteins to '
                                'species specified in the species '
                                'mapping file.')

        mandatory_arguments = parser.add_argument_group(title=''
                                                        'Required arguments')
        input_group = mandatory_arguments.add_mutually_exclusive_group(
            required=True)
        input_group.add_argument('-i', '--id',
                                 type=lambda p: Path(p).resolve(),
                                 help='Path to a text file containing '
                                 'protein IDs')
        input_group.add_argument('-f', '--fasta',
                                 type=lambda p: Path(p).resolve(),
                                 help='Path to a text file containing '
                                 'proteins in '
                                 'FASTA format')

        mandatory_arguments.add_argument('-c', '--config',
                                         type=lambda p: Path(p).resolve(),
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

    # The config argument is required
    config_file = arguments.config
    protein_params = set_params(config_file)

    """ Load the mapping file which translates species IDs between databases.
    """
    spec_mappings = species_mapping(protein_params)

    only_update_all_distances = False
    if arguments.distance:
        only_update_all_distances = True

    # This is a special species ID, where no protein ID is needed and its only
    # purpose is to update the distances between all species in the species
    # list. Prottrace is exited when the distance calculation is finished.
    if only_update_all_distances:
        protein_params.species = 'ALL'
    calculate_species_distances(protein_params, spec_mappings)

    # One of these arguments is required already. Both gen functions generate
    # the same basic information. Fasta arguments just also fill in the query
    # sequence.
    if arguments.id:
        # Adds the query sequences via the config.
        query_list = gen_ids(arguments.id, spec_mappings, protein_params)
        if query_list is None:
            print_error('An error occurred when reading query IDs')
    if arguments.fasta:
        query_list = gen_fasta(arguments.fasta, spec_mappings)
        if query_list is None:
            print_error('An error occurred when reading query IDs')

    process_query_list(query_list, protein_params)


def gen_ids(id_file, spec_mapping, config):
    """ Reads the id file into protein ID objects. """

    # You can only specify one species within the config as your query species.
    # Therefore, we assume that all query proteins represent this species and
    # that we can load its proteome once to fill in the query sequences.
    query_proteome = proteome(config, spec_mapping.get_species(config.species))

    try:
        with id_file.open('r') as i_file:
            for line in i_file:
                query_id = line.rstrip()
                yield prot_id(query_id, spec_mapping,
                              query_proteome.get_protein(query_id))
    except KeyError:
        print_error(f'The query species "{config.species}" is not the species '
                    f'of the given query protein ID "{query_id}"')


def gen_fasta(fasta_file, spec_mapping):
    """ Reads the fasta_file into records. """

    for record in fasta.gen_fasta_records(fasta_file):
        yield prot_id(record.id, spec_mapping, record.seq)


def process_query_list(id_list, protein_params):
    """ Runs all requested ProtTrace steps for a given id_list. """
    for query in id_list:
        print_progress(f'Running for ID: {query.id}')
        if protein_params.preprocessing:
            evol_model(query, protein_params)
        # Trace the detectability of the query protein.
        if protein_params.traceability_calculation:
            trace(query, protein_params)
        # Calculate the evolutionary traceability from raw trace values.
        calc_trace(query, protein_params)
        if protein_params.colorize_species_tree:
            colorize_tree(query, protein_params)


if __name__ == "__main__":
    time_report.format_verbose('Start time')
    time_record = time_report()

    # Execute the software.
    main()

    time_record.print_time('TOTAL TIME', 'hrs')
    time_report.format_verbose('End time')
