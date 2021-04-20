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

import sys
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from databases.oma_api.api import oma_sq, oma_pw, oma_gr
from databases.fdog_api.api import fdog_api, fdog_api_species

from utils.configure import set_params
from utils.file import generate_splitted_lines
from utils.log import print_error


class prot_id:
    __slots__ = ['id', 'spec_names', 'seq']

    def __init__(self, pr_id, spec_mapping, seq=None):
        self.id = pr_id
        self.spec_names = spec_mapping
        self.seq = seq

    def spec_name(self):
        return self.spec_names.species_id

    def spec_equals(self, other, form=None):
        """ Returns whether the other species name is equal to this protein's
        species name. """

        if form is None:
            return self.spec_names.species_id == other
        elif form == 'oma':
            return self.spec_names.oma == other
        elif form == 'ncbi':
            return self.spec_names.ncbi == other

    def filepath(self):
        """ Generates a fasta filename that remains valid as long as the
        protein id remains the same. """
        return Path(f'seq_{self.id}.fa')

    def write_to_fasta(self):
        """ Writes a fasta file containing its protein ID and its sequence.
        Returns the name of the written file. """
        filename = self.filepath()
        if not filename.exists():
            record = SeqRecord(self.seq, id=self.id)
            SeqIO.write(record, filename, 'fasta')
        return filename

    def remove_fasta(self):
        """ Deletes a fasta file made from its own protein ID. """
        fasta_path = Path(self.filepath())
        fasta_path.unlink(missing_ok=True)


class species_mapping_unit:
    __slots__ = ['name', 'species_id', 'fdog', 'ncbi', 'oma']

    def __init__(self, mapping_columns):
        self.name = mapping_columns[1]
        if self.column_exists(mapping_columns, 0):
            self.fdog = mapping_columns[0]
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
        with config.species_map.open('r') as sm:
            for line in generate_splitted_lines(sm):
                yield species_mapping_unit(line)

    def __str__(self):
        return f'Name: {self.name}\tID: {self.species_id}\tNCBI: {self.ncbi}\t'
        'OMA: {self.oma}'


class species_mapping:
    __slots__ = ['mappings']

    def __init__(self, config):
        self.mappings = list(species_mapping_unit.
                             generate_species_mapping(config))

    def all_species(self):
        return [unit.id for unit in self.mappings]

    def combinations(self):
        for unit_1 in self.mappings:
            for unit_2 in self.mappings:
                if unit_1 != unit_2:
                    yield (unit_1, unit_2)

    def get_species(self, species):
        for id_mapping in self.mappings:
            if id_mapping.species_id == species:
                return id_mapping

    def get_spec_from_prot(self, protein):
        oma_spec = oma_sq.prot_id_to_spec_id(protein)
        return self.get_species(oma_spec)

    def __str__(self):
        return '\n'.join([str(m) for m in self.mappings])


class orth_group:
    __slots__ = ['members', 'path']

    def __init__(self, config, prot_name, proteomes):
        mapping = species_mapping(config)
        group_api = self.resolve_group_api(config)
        self.members = list(self.collect_members(prot_name,
                                                 group_api,
                                                 mapping))
        self.collect_sequences(proteomes, config)
        self.path = self.filepath(prot_name)

    def collect_members(self, prot_name, api, mapping):
        for member_prot in api.gen_members(prot_name):
            yield prot_id(member_prot, mapping.get_spec_from_prot(member_prot))

    def collect_sequences(self, preloaded_proteomes, config):
        """ Get the sequences of orthologous group members. Previously
        instantiated proteomes are searched for first, before loading other
        species' proteomes from fasta files. """
        for member in self.members:
            for p in preloaded_proteomes:
                if p.species_name == member.spec_name():
                    member.seq = p.get_protein(member.id)
            # If no preloaded proteome is available, load others on-demand.
            proteome_on_demand = proteome(config, member.spec_name())
            member.seq = proteome_on_demand.get_protein(member.spec_name())
            # The just loaded proteome is added to the temporary list of
            # preloaded proteomes. They are lost when exiting this function.
            preloaded_proteomes.append(proteome_on_demand)

    def at_least_4_sequences(self):
        return len(self.members) > 3

    def filepath(self, prot_name):
        """ Generates a valid filename as long as the protein name stays
        the same. """
        return Path(f'ogSeqs_{prot_name}.fa')

    def write_to_file(self):
        with self.path.open('w') as orth_file:
            for member in self.members:
                member_record = SeqRecord(
                    member.seq,
                    id=member.id)
                SeqIO.write(member_record, orth_file, 'fasta')

    def resolve_group_api(self, config):
        if config.search_oma_database:
            return oma_gr
        else:
            return fdog_api


class fasta:
    """ A container of SeqIO records. """
    __slots__ = ['records']

    def __init__(self, fasta):
        self.records = list(SeqIO.parse(str(fasta)))

    @staticmethod
    def gen_fasta_records(filename):
        for record in SeqIO.parse(str(filename)):
            yield record


class proteome:
    __slots__ = ['source_api', 'species_name']

    def __init__(self, config, species_name):
        if species_name == '':
            print_error('Access to a genome stored in HaMStR requires a '
                        'species name!')
            sys.exit()
        self.species_name = species_name
        self.source_api = self.resolve_api(config, species_name)
        if not self.source_api.seqs_file.exists():
            print_error(f'The fasta file {self.source_api.seqs_file} is '
                        'invalid')

    def resolve_api(self, config, species_name):
        """ Resolves the preferred source of proteomes. The respective api
        class is initialized and returned. """
        if config.search_oma_database:
            return oma_sq(config, species_name)
        else:
            return fdog_api_species(config, species_name)

    def get_protein(self, identifier):
        """ Returns the sequence record of the given identifier. """
        return self.source_api.get_protein(identifier)

    # def get_index(self):
    #     """ Loads the proteome sequence dictionary and returns it. It is
    #     called an index, because it is typically loaded with SeqIO.index. """
    #     self.source_api.load_index()
    #     return self.source_api.index


def gen_proteomes(config, *species):
    """ Generates a proteome api for every given species. """

    for s in species:
        yield proteome(config, s)


def gen_pairwise_orthologs(config, species_1, species_2):
    """ Accesses more specific APIs to get pairwise orthologs between species.
    """

    pw_api = oma_pw(config)
    return pw_api.get_pairwise_orthologs(species_1, species_2)


if __name__ == "__main__":
    """ This script is only executed as a script for testing purposes! """
    config_file = Path('test.config')
    config = set_params(config_file)

    test_mapping = species_mapping(config)
    print(test_mapping)
