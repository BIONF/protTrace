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

""" This module manages access to OMA database information. """

import sys
from pathlib import Path

from Bio.SeqIO import index as SeqIO_index

from utils.file import generate_splitted_lines
from utils.log import print_warning, print_error


class oma:
    """ Manages content from OMA database. """

    @staticmethod
    def prot_id_to_spec_id(protein_id):
        """ Extracts the OMA species ID from an OMA protein ID. """
        return protein_id[0:5]

    @staticmethod
    def prot_ids_to_spec_ids(protein_ids):
        """ Extracts the OMA species IDs from OMA protein IDs. """
        for prot in protein_ids:
            yield prot[0:5]


class oma_pw(oma):
    """ Manages OMA pairwise orthologs. """
    __slots__ = ['oma_pairs']

    def __init__(self, config):
        self.oma_pairs = config.path_oma_pairs

    def get_pairwise_orthologs(self, species_1, species_2):
        """ Lazily reads the lines of pairwise ortholog protein IDs of the
        appropriate species pair. """

        def generate_pairs(oma_pairs, species_1, species_2):
            """ Collects pairwise orthologous protein ID lines from the
            species-specific pairwise protein ID table. """

            # Species can be stored as any name combination, separated
            # with underscores.
            sp_oma_pairs_1 = Path(f'{oma_pairs}/{species_1}_{species_2}.txt')
            sp_oma_pairs_2 = Path(f'{oma_pairs}/{species_2}_{species_1}.txt')

            if sp_oma_pairs_1.exists():
                definitive_path = sp_oma_pairs_1
            elif sp_oma_pairs_2.exists():
                definitive_path = sp_oma_pairs_2
            with definitive_path.open('r') as pair_src:
                for line in pair_src:
                    yield line

        return generate_pairs(self.oma_pairs, species_1, species_2)


class oma_gr(oma):
    __slots__ = ['oma_groups']
    """ Manages OMA groups. """

    def __init__(self, config):
        self.oma_groups = config.path_oma_group

    def gen_members(self, prot_name):

        def gen_gr_lines(instance):
            with instance.oma_groups.open('r') as gr_file:
                for line in generate_splitted_lines(gr_file):
                    yield line

        for line in gen_gr_lines(self):
            if prot_name in line:
                # The first two elements of an oma_groups.txt file shows the
                # group number and fingerprint. The member protein ids are
                # shown from element 3 onwards.
                for member in line[2:]:
                    yield member


class oma_sq(oma):
    __slots__ = ['seqs_file', 'index']
    """ Manages OMA sequences. """

    def __init__(self, config, species_name):
        if not config.path_oma_seqs.is_dir():
            print_error('The oma seqs path in the config must lead to a '
                        'directory!')
        self.seqs_file = (config.path_oma_seqs /
                          species_name).with_suffix('.fa')
        self.index = self.index_seqs()

    def get_protein(self, identifier):
        """ Searches the protein sequence in the SeqRecord index of the
        oma-seqs file. """
        return str(self.index[identifier].seq)

    def get_proteome_file(self):
        """ Returns the Path of the original proteome that this class has
        been instantiated with. The proteome must be exclusive to one species.
        """
        return self.seqs_file

    def index_seqs(self, necessary_ids=None):
        try:
            sequence_index = SeqIO_index(str(self.seqs_file), 'fasta')
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
