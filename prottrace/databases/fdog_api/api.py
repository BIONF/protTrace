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

""" This module manages access to fDog interfaces. """

from pathlib import Path
from Bio import SeqIO


class fdog_api():
    """ Manages fdog related functions. """
    __slots__ = ['path', 'environment']

    def __init__(self, prot_config):
        """ Class initialiser """
        # The path to the fdog main directory.
        self.path_main = prot_config.fdog
        self.path_genome = self.path_main / 'fdog' / 'data' / 'genome_dir'

    def get_genome(self, species_name):
        return Path(self.path_genome / species_name)


class fdog_api_species(fdog_api):
    """ Accesses the genome of a particular species in fDOG. """
    __slots__ = ['species_name']

    def __init__(self, prot_config, species_name):
        super().__init__(prot_config)
        self.species_name = species_name

    def get_protein(self, identifier):
        prot = super(fdog_api_species, self).get_genome(self.species_name)
        for record in SeqIO.parse(str(prot)):
            if record.id == identifier:
                return record
