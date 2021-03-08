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

# Module to read in program config file and prepare variables
# to be imported in other modules

import os
import sys
from prottrace.utils.log import print_error


class setParams:
    """ Stores the configuration for ProtTrace. """

    def string_option(self, value):
        """ Parses the value of the name value pair
            to the searched name. """
        return str(value)

    def integer_option(self, value):
        """ Parses the value of the name value pair
            to the searched name into an integer. """
        return int(value)

    def float_option(self, value):
        """ Parses the value of the name value pair
            to the searched name into a float. """
        return float(value)

    def boolean_option(self, value):
        """ Checks the name value pair for YES and NO and
            returns True and False, depending on the
            result. """
        if str(value) in set(['YES', 'yes', 'True', 'true']):
            return True
        else:
            return False

    def path_option(self, value):
        """ Returns the absolute path to the path given
            in the name value pair. """
        return os.path.abspath(value)

    def create_path_option(self, value):
        """ Creates the path specified in the name value pair
            if it does not exist yet. """
        created_path = self.path_option(value)
        if not os.path.exists(created_path):
            os.mkdir(created_path)
        return created_path

    def read_config_file(self, config_file):
        """ Reads the config file in, splits it into lines, ignores empty and
            commented lines. """
        valid_lines = []
        with open(config_file, 'r') as config:
            for line in config:
                line = line.rstrip()
                # Rejects empty lines and possibly indented lines starting
                # with a hashtag for comments.
                if not line == '' and line.lstrip()[0] != '#':
                    valid_lines.append(line)
        return valid_lines

    def parse_key_value_pairs(self, config_file):
        """ Extracts key:value pairs out of configuration lines and checks
            for duplicates. Other general validation processes happen
            here. """

        key_value_pairs = {}
        for line in self.read_config_file(config_file):
            splitted_line = line.split(':')

            # Testing for key:value format.
            if len(splitted_line) != 2:
                raise ValueError('The line \n{0}\n does not contain exactly '
                                 'one pair of type key:value!'
                                 .format(splitted_line))

            # Testing for duplicates.
            if splitted_line[0] in key_value_pairs:
                raise ValueError('The configuration key {0} is duplicated!'
                                 .format(splitted_line[0]))

            key_value_pairs[splitted_line[0]] = splitted_line[1]
        return key_value_pairs

    def __init__(self, config_file):

        configs = self.parse_key_value_pairs(config_file)
        try:
            # This option sets the query species
            self.species = self.string_option(configs['species'])
            # Yes-No settings
            self.search_oma_database = self.boolean_option(
                configs['search_oma_database'])
            self.orthologs_prediction = self.boolean_option(
                configs['orthologs_prediction'])
            self.run_hamstr = self.boolean_option(configs['run_hamstr'])
            self.run_hamstrOneSeq = self.boolean_option(
                configs['run_hamstrOneSeq'])
            if not self.run_hamstrOneSeq:
                self.run_hamstr = False
            self.fas_score = self.boolean_option(configs['fas_score'])
            self.preprocessing = self.boolean_option(configs['preprocessing'])
            self.traceability_calculation = self.boolean_option(
                configs['traceability_calculation'])
            self.calculate_scaling_factor = self.boolean_option(
                configs['calculate_scaling_factor'])
            self.calculate_indel = self.boolean_option(
                configs['calculate_indel'])
            self.perform_msa = self.boolean_option(configs['perform_msa'])
            self.delete_temp = self.boolean_option(
                configs['delete_temporary_files'])
            self.reuse_cache = self.boolean_option(configs['reuse_cache'])
            self.mapTraceabilitySpeciesTree = self.boolean_option(
                configs['map_traceability_tree'])
            self.includeParalogs = self.boolean_option(
                configs['include_paralogs'])
            self.phylogeneticTreeReconstruction = self.boolean_option(
                configs['orthologs_tree_reconstruction'])
            self.run_spartaABC = self.boolean_option(
                configs['run_spartaABC'])
            self.evolve_dawg = self.boolean_option(
                configs['dawg_instead_of_indelible'])

            # Options where numbers or strings can be specified
            self.aa_substitution_matrix = self.string_option(
                configs['aa_substitution_matrix'])
            self.simulation_runs = self.integer_option(
                configs['simulation_runs'])
            self.nr_processors = self.integer_option(
                configs['nr_of_processors'])

            # Default evolution simulation parameters
            self.default_indel = self.float_option(configs['default_indel'])
            self.default_indel_distribution = self.float_option(
                configs['default_indel_distribution'])
            self.default_scaling_factor = self.float_option(
                configs['default_scaling_factor'])

            # Software paths
            self.sparta = self.path_option(configs['spartaABC'])
            self.msa = self.path_option(configs['linsi'])
            self.REvolver = self.path_option(configs['REvolver'])
            self.hmmfetch = self.path_option(configs['hmmfetch'])
            self.hmmscan = self.path_option(configs['hmmscan'])
            self.iqtree = self.path_option(configs['iqtree'])
            self.treepuzzle = self.path_option(configs['treepuzzle'])
            # self.clustalw = self.path_option(configs['clustalw'])
            self.blastp = self.path_option(configs['blastp'])
            self.makeblastdb = self.path_option(configs['makeblastdb'])
            self.R = self.path_option(configs['Rscript'])

            # HaMStR
            self.hamstr = self.path_option(configs['hamstr'])
            self.hamstrOneSeq = self.path_option(configs['oneseq'])
            self.hamstr_environment = self.string_option(
                configs['hamstr_environment'])
            if self.hamstr_environment == 'default':
                self.hamstr_environment = ''

            # Set output and cache directory paths or create default ones
            # Each protein creates its own subdirectory within this output
            # directory
            self.path_work_dir = self.create_path_option(
                configs['path_output_dir'])
            # The cache directory contains global intermediary files,
            # such as calculated species distances
            self.path_cache = self.create_path_option(configs['path_cache'])
            # Pairwise species distances are calculated within this directory
            self.path_distance_work_dir = self.create_path_option(
                configs['path_distances'])

            # Input file paths
            self.species_MaxLikMatrix = self.path_option(
                configs['species_MaxLikMatrix'])
            self.hamstr_oma_tree_map = self.path_option(
                configs['Xref_mapping_file'])
            self.species_MaxLikMatrix = self.path_option(
                configs['species_MaxLikMatrix'])
            self.path_oma_seqs = self.path_option(configs['path_oma_seqs'])
            self.path_oma_pairs = self.path_option(configs['path_oma_pairs'])
            self.path_oma_group = self.path_option(configs['path_oma_group'])
            self.pfam_database = self.path_option(configs['pfam_database'])
            self.species_tree = self.path_option(
                configs['reference_species_tree'])
            self.simulation_tree = self.path_option(configs['simulation_tree'])
            self.concat_alignments_script = self.path_option(
                configs['concat_alignments'])
            self.decay_script = self.path_option(configs['decay_script'])
            self.plot_figtree = self.path_option(configs['plot_figtree'])
            self.fas_annotations = self.path_option(configs['fas_annotations'])
        except KeyError as e:
            print_error('Configuration key "{0}" not found!'.format(e.args[0]))
            sys.exit()
