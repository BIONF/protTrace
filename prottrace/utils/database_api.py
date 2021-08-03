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

""" This module offers abstract classes to implement new APIs for access to
the data of databases or orthology prediction tools. These functions should be
seen as guidelines to predict the calls of the data_api module.  """


class database_pw:
    """ Provides an interface to offer pairwise orthologs to classes of the
    data_api module. """
    __slots__ = []

    def get_pairwise_orthologs(self, species_1, species_2):
        """ Returns an iterable tuple of two pairwise orthologous protein
        IDs. """
        raise NotImplementedError('This database does not offer pairwise '
                                  'orthologs!')


class ortholog_info:
    """ Since orthologous group retrieval methods might vary, information to
    the data_api module should be returned as instances of this container
    class. """
    __slots__ = ['member_id', 'member_seq', 'species_id', 'database']

    def __init__(self, member_id, member_seq=None, species_id=None,
                 database=None):
        # The orthologous group member protein ID MUST be provided
        self.member_id = member_id

        # The following values ease the processing by data_api, but may not
        # be necessary.
        self.member_seq = member_seq
        self.species_id = species_id
        self.database = database


class database_gr:
    """ Provides an interface to offer orthologous groups to classes of the
    data_api module. """
    __slots__ = []

    def gen_members(self, query):
        """ Generates an iterable of ortholog_info instances. query is an
        instance of data_api.prot_id. """
        raise NotImplementedError('This database does not offer orthologous '
                                  'groups!')


class database_sq:
    """ Provides an interface to offer protein sequences to classes of the
    data_api module. """
    __slots__ = []

    def get_protein(self, identifier):
        """ Lookup the FASTA file or the indexed FASTA file for the requested
        identifer. """
        raise NotImplementedError('This database does not offer the sequence '
                                  'of single proteins!')

    def get_proteome_file(self):
        """ Get the FASTA file, where this database stores its sequences. """
        raise NotImplementedError('This database does not offer full FASTA '
                                  'files!')

    def index_seqs(self):
        """ Uses Bio.SeqIO.index to map the proteome FASTA file into a
        dictionary. """
        raise NotImplementedError('This database does not offer indexed '
                                  'FASTA file access!')
