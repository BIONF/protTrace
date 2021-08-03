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
from subprocess import (
    run,
    CalledProcessError
)
from Bio.SeqIO import parse as SeqIO_parse

from utils.database_api import (
    ortholog_info,
    database_gr
)
from utils.log import print_error, print_progress


class fdog():
    """ Manages fdog related functions. """
    __slots__ = ['path', 'environment']

    def __init__(self, config):
        """ Class initialiser """
        # The path to the fdog main directory.
        self.path_main = config.fdog
        self.fdog_envir = config.fdog_environment
        self.path_environ = self.resolve_environment(self.fdog_envir)
        self.path_genome = (self.path_main / 'fdog' / 'data'
                            / 'genome_dir' + self.path_environ)

    def get_genome(self, species_name):
        return Path(self.path_genome / species_name)

    def resolve_environment(envir):
        """ Prepends an underscore to the fdog environment variable, if it
        has been specified. """

        if envir == '':
            return envir
        else:
            return '_' + envir


class fdog_gr(fdog, database_gr):
    """ Manages Fdog orthologous groups. """
    __slots__ = ['del_temp', 'include_paralogs', 'cores']

    def __init__(self, config):
        fdog.__init__(config)
        self.del_temp = config.delete_temp
        self.include_paralogs = config.include_paralogs
        self.cores = config.nr_processors

    def gen_members(self, query):
        """ Data API interface to retrieve orthologous group member IDs. """

        self.one_seq(self, query)

        return self.record_fdog_output(query)

    def one_seq(self, query):
        """ Assumes that the api between HaMStR One-Seq and Fdog has not
        changed. """

        # Check whether the query species is actually available in FDOG.
        if not self.get_genome(query.spec_names.fdog).exists():
            print_error('The query species does not exist among the fDOG '
                        'genomes')

        # fDOG needs the fasta file of the query protein. We write the query
        # sequence into a temporary file.
        ref_seq_file = query.write_to_fasta()

        fdog_command = ['fdog.run',
                        f'-seqFile={str(ref_seq_file)}',
                        f'-seqName={str(query.id)}',
                        f'-refspec={str(query.spec_names.fdog)}',
                        f'-cpu={str(self.cores)}',
                        '-coreOrth=5',
                        '-minDist=genus',
                        '-maxDist=superkingdom',
                        '-checkCoorthologsRef',
                        '-fasoff',
                        '-local',
                        '-strict']

        # Variable arguments to fDOG/HaMStR.
        if self.fdog_envir != '':
            fdog_command.append(f'-addenv={str(self.fdog_envir)}')
        if self.del_temp:
            fdog_command.append('-cleanup')
        if Path('inputTaxaSet_oneSeq.txt').exists():
            fdog_command.append('-coreTaxa=inputTaxaSet_oneSeq.txt')
        if self.include_paralogs:
            fdog_command.append('-rep')

        print_progress('Running fDOG.')

        try:
            run(fdog_command, check=True)
        except CalledProcessError as e:
            print_error('fDOG exited with an error')
            raise e

    def record_fdog_output(self, query):
        """ Notices the orthologous group members in the extended core
        orthologs file from fDOG and returns the member's IDs. """

        fdog_output = Path(query.id + '.extended.fa')
        with fdog_output.open('r') as extended:
            for line in extended:
                if '>' in line:
                    splitted_line = line.split('|')
                    fdog_ID = splitted_line[1]
                    seq_identifier = splitted_line[2]
                    seq = (next(extended).rstrip().
                           replace('*', '').replace('X', ''))
                    yield ortholog_info(seq_identifier, seq, fdog_ID, 'fdog')


class fdog_api_species(fdog):
    """ Accesses the genome of a particular species in fDOG. """
    __slots__ = ['species_name']

    def __init__(self, config, species_name):
        fdog.__init__(config)
        self.species_name = species_name

    def get_protein(self, identifier):
        prot = super(fdog_api_species, self).get_genome(self.species_name)
        for record in SeqIO_parse(str(prot)):
            if record.id == identifier:
                return record
