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

""" This python module organizes data access for ProtTrace. """

import os
import sys
import argparse
import prottrace.utils.configure
from collections import OrderedDict
from pathlib import Path
from prottrace.utils.log import print_warning, print_error
from Bio import SeqIO


class species_mapping_unit:
    __slots__ = ['name', 'species_id', 'ncbi', 'oma']

    def __init__(self, mapping_columns):
        self.name = mapping_columns[1]
        if self.column_exists(mapping_columns, 0):
            self.hamstr_id = mapping_columns[0]
        if self.column_exists(mapping_columns, 2):
            self.ncbi = mapping_columns[2]
            self.species_id = mapping_columns[2]
        if self.column_exists(mapping_columns, 3):
            self.oma = mapping_columns[3]
            self.species_id = mapping_columns[3]

    def column_exists(self, mapping_columns, index):
        return (len(mapping_columns) >= index
                and mapping_columns[index] != '')

    @staticmethod
    def generate_species_mapping(config):
        with open(config.hamstr_oma_tree_map, 'r') as sm:
            for line in sm:
                yield species_mapping(line.rstrip().split('\t'))


class species_mapping:
    __slots__ = ['mappings']

    def __init__(self, config):
        self.mappings = list(species_mapping_unit.
                             generate_species_mapping(config))

    def get_species(self, species):
        for id_mapping in self.mappings:
            if id_mapping.species_id == species:
                return id_mapping


class proteome:
    __slots__ = ['source_name', 'source_file', 'species_name']

    def __init__(self, config, species_name=''):
        self.species_name = species_name
        if config.search_oma_database:
            self.source_file = config.path_oma_seqs
        else:
            if self.species_name == '':
                print_error('Access to a genome stored in HaMStR requires a '
                            'species name!')
                sys.exit()
            self.source_file = self.get_hamstr_gen_path(config,
                                                        self.species_name)
        # The file basename is stored for sanity checks.
        self.source_name = self.source_file.split('/')[-1]

    def get_hamstr_gen_path(self, config, species_name):
        """ Builds the path to a genome stored within the HaMStR directory. """
        hamstr_config = hamstr_api(config)
        return hamstr_config.resolve_path('genome_dir') + '/' + species_name

    def get_indexed_proteome(self):
        """ SeqIO can index a fasta file for dictionary-like access. This is
        useful for random access when the protein ids are known. """
        return SeqIO.index(self.source_file, 'fasta')


class hamstr_api:
    """ Manages hamstr related functions. """
    __slots__ = ['path', 'environment']

    def __init__(self, prot_config):
        """ Class initialiser """
        self.path = prot_config.hamstr
        self.environment = prot_config.hamstr_environment

    def resolve_path(self, directory):
        """ Evaluates whether the directory must be adjusted to a
        non-standard environment. """
        absolute_directory = self.path + "/" + directory
        if self.environment != "":
            absolute_directory += "_" + self.environment
            return absolute_directory


class oma_api:
    """ Manages content from OMA database. """
    __slots__ = ['oma_pairs', 'oma_seqs', 'oma_group', 'oma_pairs_seek']

    def __init__(self, config):
        self.oma_pairs = config.path_oma_pairs

        ipf = Path(str(self.oma_pairs))
        self.oma_pairs_seek = ipf.with_suffix(ipf.suffix + '.indexed')

        self.oma_seqs = config.path_oma_seqs
        self.oma_group = config.path_oma_group

    @staticmethod
    def prot_id_to_spec_id(protein_id):
        """ Extracts the OMA species ID from an OMA protein ID. """
        return protein_id[1:5]

    @staticmethod
    def prot_ids_to_spec_ids(protein_ids):
        """ Extracts the OMA species IDs from OMA protein IDs. """
        for prot in protein_ids:
            yield prot[1:5]

    def index_pairs(self):
        """ Adds seek positions to a previously sorted oma_pairs.txt file.
            The species IDs should have also been added as 3th and 4th columns.
            """
        pair_index = []

        def format_pair_index(species, pos):
            """ Produces the text output of a species pair - position dataset.
            """
            pair_format = ['%']
            pair_format.extend(species)
            pair_format.extend([str(p) for p in pos])
            return '\t'.join(pair_format)

        def pair_seek_positions(filename):
            """ Annotates the byte positions of species pairs in an
                orthologous pair file, such as oma_pairs.txt. """
            # The same set of species pair can appear at very different
            # sections of the file.
            pair_index = OrderedDict()
            previous = ''

            with open(filename, 'r') as pairs:
                # The position is told before reading the line.
                for pos, spl_line in generate_splitted_lines_with_pos(pairs):
                    # frozenset enables their use as dict keys
                    current_species = frozenset(spl_line[2:4])
                    if current_species != previous:
                        if previous != '':
                            # The previous ending position tells us the length
                            # of the section
                            pair_index[previous].append(
                                pos - pair_index[previous][-1])

                        # The starting position of the current species pair is
                        # added.
                        if current_species not in pair_index:
                            # OrderedDict preserves the order.
                            pair_index[current_species] = []
                        pair_index[current_species].append(pos)

                        # Update the previous species for the upcoming section.
                        previous = current_species

            return '\n'.join([format_pair_index(key, pair_index[key])
                              for key in pair_index]) + '\n'

        pair_index = pair_seek_positions(self.oma_pairs)

        with self.oma_pairs_seek.open('w', encoding='utf-8') as indxd_pairs:
            for ind in pair_index:
                indxd_pairs.write(ind)

    def get_pairwise_orthologs(self, species_1, species_2):
        """ Lazily reads the lines of pairwise ortholog protein IDs of the
        appropriate species pair. """

        def search_seek_pnts_legacy(filepath, checked_species):
            """ Collect the seek points for species of interest. """
            seek_points = []
            with filepath.open('r') as handler:
                for spl_line in generate_splitted_lines(handler):
                    if spl_line[0] == '%':
                        if set(spl_line[1:3]) == checked_species:
                            seek_points.append([int(pnt) for pnt
                                                in spl_line[3:]])
                    else:
                        break
            return seek_points

        def generate_pairs_legacy(filename, seek_points, species_1, species_2):
            with open(filename, 'r', encoding='utf-8') as handler:
                for pnts in seek_points:
                    for i in range(0, len(pnts), 2):
                        handler.seek(pnts[i], os.SEEK_SET)
                        content = handler.read(pnts[i+1])
                        content = content.split('\n')
                        # The last bit can be empty because the pointer shows
                        # where to exactly start the next section.
                        if content[-1] == '':
                            del content[-1]

                        # Validation of the first entry.
                        first_entry = content[0].split('\t')
                        if (not set(first_entry[2:4]) ==
                                set((species_1, species_2))
                                and len(first_entry) == 4):
                            print_error('The pairwise orthologous section was '
                                        'wrong! Species pair: {0} {1}'
                                        .format(species_1, species_2))
                            raise ValueError

                        # Validation of the last entry.
                        last_entry = content[-1].split('\t')
                        if (not set(last_entry[2:4]) ==
                                set((species_1, species_2))
                                and len(last_entry) == 4):
                            print_error('The pairwise orthologous section was '
                                        'wrong! Species pair: {0} {1}'
                                        .format(species_1, species_2))
                            raise ValueError
                            print_error('The pairwise orthologous section was '
                                        'wrong! Species pair: {0} {1}'
                                        .format(species_1, species_2))
                            raise ValueError

                        # The content is given out per line.
                        for line in content:
                            yield line

        def generate_pairs(oma_pairs, species_1, species_2):
            """ Collects pairwise orthologous protein ID lines from the
            species-specific pairwise protein ID table. """

            sp_oma_pairs_1 = '{0}/{1}_{2}.txt'.format(oma_pairs, species_1,
                                                      species_2)

            sp_oma_pairs_2 = '{0}/{2}_{1}.txt'.format(oma_pairs, species_1,
                                                      species_2)
            if os.path.exists(sp_oma_pairs_1):
                definitive_path = sp_oma_pairs_1
            elif os.path.exists(sp_oma_pairs_2):
                definitive_path = sp_oma_pairs_2
            with open(definitive_path, 'r') as pair_src:
                for line in pair_src:
                    yield line

        # query_species = set((species_1, species_2))
        # seek_points = search_seek_pnts(self.oma_pairs_seek, query_species)
        return generate_pairs(self.oma_pairs, species_1, species_2)

    def index_seqs(self, necessary_ids=None):
        try:
            sequence_index = SeqIO.index(self.oma_seqs, 'fasta')
            # The produced index is checked for all necessary species.
            if necessary_ids is not None:
                self.index_seqs_check(sequence_index, necessary_ids)
            return sequence_index
        except KeyboardInterrupt:
            print_warning('User interrupted the indexing of OMA sequences.')
            sys.exit()
        except FileNotFoundError as e:
            print_error('The OMA seqs FASTA file is missing!')
            raise e

    def index_seqs_check(self, sequence_index, identifiers):
        """ Returns an indexed oma_seqs file after checking whether all
        provided OMA species IDs are present. """

        sequence_index = self.index_seqs()
        for ident in identifiers:
            if ident not in sequence_index:
                print_error('{0} could not be found in the OMA seqs file!')
                sys.exit()


def index_oma_pairs(config, force=False, testpair=None):
    config = configure.setParams(arguments.config)

    oma = oma_api(config)
    if not os.path.exists(oma.oma_pairs_seek) or force:
        oma.index_pairs()

    if testpair is not None:
        test_spec = arguments.testpair.split('_')
        return oma.get_pairwise_orthologs(test_spec[0], test_spec[1])
    return 0


def generate_splitted_lines_with_pos(handler):
    while True:
        previous_pos = handler.tell()
        line = handler.readline()
        if not line:
            break
        yield (previous_pos, line.rstrip().split('\t'))


def generate_splitted_lines(lines):
    while True:
        line = lines.readline()
        if not line:
            break
        yield line.rstrip().split('\t')


if __name__ == "__main__":
    def Argparse():
        parser = argparse.ArgumentParser(description='Prepares the larger '
                                         'sized OMA pairs file for reading in '
                                         'pairwise orthologs.')
        parser.add_argument('-c', '--config', type=str, help='The path to a '
                            'ProtTrace configuration file.')
        parser.add_argument('-f', '--force', action='store_true', help='Force '
                            'overrides an existing index file.')
        parser.add_argument('-t', '--testpair', type=str, nargs='?',
                            default=None, help='Define a '
                            'pair of "species1_species2" to get from the '
                            'indexed OMA pairs file.')
        return parser.parse_args()

    arguments = Argparse()

    outcode = 1
    if arguments.testpair is not None:
        print(index_oma_pairs(arguments.config, arguments.force,
                              arguments.testpair))
        outcode = 0
    else:
        outcode = index_oma_pairs(arguments.config, arguments.force,
                                  arguments.testpair)

    sys.exit(outcode)
