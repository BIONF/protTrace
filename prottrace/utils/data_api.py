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

from Bio.SeqIO import (
    parse as SeqIO_parse,
    write as SeqIO_write
)
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from databases.oma_api.api import oma_sq, oma_pw, oma_gr
from databases.fdog_api.api import fdog_api, fdog_api_species

from utils.configure import set_params
from utils.file import generate_splitted_lines as generate_lines
from utils.log import (
    print_warning,
    print_error
)


class prot_id:
    __slots__ = ['id', 'spec_names', 'seq']

    def __init__(self, pr_id, spec_mapping, seq=None):
        self.id = pr_id
        if isinstance(spec_mapping, species_mapping):
            self.spec_names = spec_mapping.get_spec_from_prot(pr_id)
        else:
            self.spec_names = spec_mapping
        self.seq = seq

    @property
    def spec_id(self):
        try:
            return self.spec_names.species_id
        except AttributeError:
            print_error('The species mapper of this protein could not be '
                        'found')
            sys.exit()

    @property
    def spec_name(self):
        try:
            return self.spec_names.taxon
        except AttributeError:
            print_error('The species name of this protein could not be found')
            sys.exit()

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
        return Path(f'seq_{self.id}.fa').resolve()

    def write_to_fasta(self):
        """ Writes a fasta file containing its protein ID and its sequence.
        Returns the name of the written file. """
        filename = self.filepath()
        if not filename.exists():
            record = SeqRecord(Seq(str(self.seq)), id=self.id,
                               description=self.spec_id)
            SeqIO_write(record, str(filename), 'fasta')
        return filename

    def remove_fasta(self):
        """ Deletes a fasta file made from its own protein ID. """
        self.filepath().unlink()

    def __str__(self):
        return f'{self.id}, {self.spec_name}'


class species_mapping_unit:
    __slots__ = ['name', 'species_id', 'fdog', 'ncbi', 'oma', 'taxon']

    def __init__(self, mapping_columns):
        # if __debug__:
        #     print('Filling mapping information for species '
        #           f'{mapping_columns[1]}')
        #     print(f'Number of mapping columns: {len(mapping_columns)}')
        self.name = mapping_columns[1]
        if self.column_exists(mapping_columns, 0):
            self.fdog = mapping_columns[0]
        if self.column_exists(mapping_columns, 1):
            self.taxon = mapping_columns[1]
        if self.column_exists(mapping_columns, 2):
            self.ncbi = mapping_columns[2]
            self.species_id = mapping_columns[2]
        if self.column_exists(mapping_columns, 3):
            self.oma = mapping_columns[3]
            self.species_id = mapping_columns[3]

    # This property provides the same id interface for species_mapping_unit as
    # for prot_id.
    @property
    def id(self):
        return self.species_id

    def match_any_id(self, query):
        return (self.match_oma(query)
                or self.match_ncbi(query)
                or self.match_fdog(query))

    def match_oma(self, query):
        return query == self.oma

    def match_ncbi(self, query):
        return query == self.ncbi

    def match_fdog(self, query):
        return query == self.fdog

    def match_id(self, query, data):
        return (data == 'oma' and self.match_oma(query)
                or data == 'ncbi' and self.match_ncbi(query)
                or data == 'fdog' and self.match_fdog(query))

    def column_exists(self, mapping_columns, index):
        return (len(mapping_columns) >= index and mapping_columns[index] != '')

    @staticmethod
    def generate_species_mapping(config):
        with config.species_map.open('r') as sm:
            for line in generate_lines(sm):
                yield species_mapping_unit(line)

    def __str__(self):
        return f'Name: {self.name}\tID: {self.species_id}\tNCBI: {self.ncbi}\t'
        'OMA: {self.oma}\n'


class species_mapping:
    """ Manages a collection of species_mapping_units. """
    __slots__ = ['mappings']

    def __init__(self, config):
        self.mappings = list(species_mapping_unit.
                             generate_species_mapping(config))

    def all_species(self):
        """ Returns a list of all species ids within this collection of
        species mappings. """
        return [unit.species_id for unit in self.mappings]

    def combinations(self):
        for unit_1 in self.mappings:
            for unit_2 in self.mappings:
                if unit_1 != unit_2:
                    yield (unit_1, unit_2)

    def get_species(self, species, database=None):

        """ Returns the species_mapping_unit of a particular species.
        Since the species mapping file dictates which species to analyze,
        we drop species that do not match here. """

        resulting_mapping = None
        for id_mapping in self.mappings:
            if database is not None:
                if database == 'oma':
                    if id_mapping.match_oma(species):
                        resulting_mapping = id_mapping
                        break
            if id_mapping.match_any_id(species):
                resulting_mapping = id_mapping
                break
        return resulting_mapping

    def get_spec_from_prot(self, protein):
        oma_spec = oma_sq.prot_id_to_spec_id(protein)
        return self.get_species(oma_spec, 'oma')

    def __str__(self):
        return '\n'.join([str(m) for m in self.mappings])


class orth_group:
    __slots__ = ['members', 'path']

    def __init__(self, config, prot_name):
        mapping = species_mapping(config)
        group_api = self.resolve_group_api(config)
        self.members = self.collect_members(prot_name, group_api, mapping)
        self.collect_sequences(config)
        self.path = self.filepath(prot_name)

        if not self.path.exists():
            self.write_to_file()

    def gen_member_combinations(self):
        """ Generates unique combinations of 2 different members.
        Returns a tuple of prot_id instances. """

        for spec1 in self.members:
            for spec2 in self.members:
                if spec1.id != spec2.id:
                    yield (spec1, spec2)

    def gen_species_ml_dists(self, config):
        """ Returns the pairwise ML species distance between each member in
        the orthologous group. """

        for combination in self.gen_member_combinations():
            spec_map1 = combination[0].spec_names
            spec_map2 = combination[1].spec_names
            dist = get_species_distance(config, spec_map1, spec_map2)
            yield species_pair_distance(spec_map1.species_id,
                                        spec_map2.species_id,
                                        dist)

    def member_prot_to_spec(self, protein):
        """ Tries to identify the query protein ID among members of this
        orthologous group and then returns the member species ID. """

        for member in self.members:
            if member.id == protein:
                return member.spec_id

        raise ValueError('The given protein ID is not a member of this '
                         'orthologous group')

    def create_prot_id(self, prot, mapping):
        """ A function to map the creation of protein entries. It mainly
        streamlines missing species mappings into None types which can be
        filtered from the final list. """
        p_mapping = mapping.get_spec_from_prot(prot)
        if p_mapping is None:
            return None
        else:
            return prot_id(prot, p_mapping)

    def collect_members(self, prot_name, api, mapping):
        """ Maps orthologous group members by their protein ID to their
        respective species ID mapping. They may access different types of
        data, so we collect the full mapping. """
        orth_maps = map(lambda member_prot: self.create_prot_id(member_prot,
                                                                mapping),
                        api.gen_members(prot_name))

        in_species_mapping = [prot for prot in orth_maps if prot is not None]

        return list(in_species_mapping)

    def collect_sequences(self, config):
        """ Get the sequences of orthologous group members. Previously
        instantiated proteomes are searched for first, before loading other
        species' proteomes from fasta files. """

        proteomes = []

        for member in self.members:
            for p in proteomes:
                if p.species_name == member.spec_id:
                    member.seq = p.get_protein(member.id)

            # If no preloaded proteome is available, load others on-demand.
            proteome_on_demand = proteome(config, member.spec_id)
            member.seq = proteome_on_demand.get_protein(member.id)
            # The just loaded proteome is added to the temporary list of
            # preloaded proteomes. They are lost when exiting this function.
            proteomes.append(proteome_on_demand)

    def at_least_4_sequences(self):
        return len(self.members) > 3

    def filepath(self, prot_name):
        """ Generates a valid filename as long as the protein name stays
        the same. """
        return Path(f'ogSeqs_{prot_name}.fa').resolve()

    def write_to_file(self):
        with self.path.open('w') as orth_file:
            for member in self.members:
                member_record = SeqRecord(
                    Seq(member.seq),
                    id=member.id,
                    description=member.spec_id)
                SeqIO_write(member_record, orth_file, 'fasta')

    def resolve_group_api(self, config):
        if config.search_oma_database:
            return oma_gr(config)
        else:
            return fdog_api(config)


class fasta:
    """ A container of SeqIO records. """
    __slots__ = ['records']

    def __init__(self, fasta):
        self.records = list(SeqIO_parse(str(fasta), 'fasta'))

    @staticmethod
    def gen_fasta_records(filename):
        for record in SeqIO_parse(str(filename), 'fasta'):
            yield record


class proteome:
    __slots__ = ['source_api', 'species_name']

    def __init__(self, config, species_name):
        if species_name == '':
            print_error('Access to a genome stored in HaMStR requires a '
                        'species name!')
            sys.exit()
        if isinstance(species_name, species_mapping_unit):
            self.species_name = species_name.species_id
        else:
            self.species_name = species_name
        self.source_api = self.resolve_api(config, self.species_name)
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

    def transfer_proteome(self, link_name):
        """ Creates a symbolic link from the associated proteome to the
        target filename. """
        source_file = self.source_api.get_proteome_file()
        link_name.symlink_to(source_file)

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


def map_distance_species(config, mappings, query=None):
    """ Returns the set of mappings of all species with valid distances in
    the ML table file. """

    def validate_distance(distance):
        """ Checks whether the number associated with the species pair is
        actually a valid species distance, e.g. a float. """

        try:
            float(distance)
            return True
        except ValueError:
            return False

    def gen_valid_ML_lines(config):
        """ Removes complexity out of the parent function. Returns
        lines splitted by tab, if the number on the third column correctly
        represents a pairwise ML distance. It also reduces the line to the
        species pair. """
        with config.ML_table.open('r') as ml_table:
            for line in generate_lines(ml_table):
                if validate_distance(line[2]):
                    yield line[0:2]

    found_species = set()
    for specs in gen_valid_ML_lines(config):
        if query is None or query == 'ALL':
            found_species.add(mappings.get_species(specs[0]).species_id)
            found_species.add(mappings.get_species(specs[1]).species_id)
        elif query in specs:
            # query_index retrieves the index of the query id in the
            # line. Possible positions are [0,1]. We want to add the
            # other species to our set.
            # If query index = 1: 1 - 1 = 0, the other species is at
            # index 0.
            # If query index = 0: 1 - 0 = 1, the other species is at
            # index 1.
            found_species.add(mappings.get_species(
                specs[1 - specs.index(query)]).species_id)

    # if __debug__:
    #     already_computed = '\t'.join(found_species)
    #     print(f'Precalculated species distances to: {already_computed}')

    return found_species


def get_table_distance(config, query, subject):
    """ Returns the pairwise species distance between the query and the
    subject. """

    species_set = frozenset([query, subject])

    with config.ML_table.open('r') as ml_table:
        for line in generate_lines(ml_table):
            if species_set == frozenset([line[0], line[1]]):
                return float(line[2])


def get_cached_distance(cache_path, query, subject):
    """ Searches the cache directory for the paired species distance. """

    possible_path_1 = Path(f'{query}_{subject}.lik')
    possible_path_2 = Path(f'{subject}_{query}.lik')
    definitive_path = None

    if possible_path_1.exists():
        definitive_path = possible_path_1
    elif possible_path_2.exists():
        definitive_path = possible_path_2

    if definitive_path is not None:
        with definitive_path.open('r') as lik_file:
            return float(lik_file.read().rstrip())
    else:
        return None


def cached_distance_first(config, query, subject):
    """ Searches the cache first for an existing pairwise species distance.
    Then searches the ML_table file. """

    cached_distance = get_cached_distance(config.path_cache, query, subject)

    if cached_distance is None:
        return get_table_distance(config, query, subject)
    else:
        return cached_distance


def get_species_distance(config, query_mapping, subject_mapping, default=None):
    """ Collects the pairwise species distance between the query and the
    subject. """

    if config.reuse_cache:
        result = cached_distance_first(config, query_mapping.species_id,
                                       subject_mapping.species_id)
    else:
        result = get_table_distance(config, query_mapping.species_id,
                                    subject_mapping.species_id)

    if result is None:
        print_warning('The pairwise distance between '
                      f'{query_mapping.species_id} and '
                      f'{subject_mapping.species_id} is invalid!')
    return result


def get_species_distances(config, query_mapping, subject_mappings,
                          default=None):
    """ Collects all pairwise species distances between the query and the
    subjects. """

    # Asks for the reuse_cache setting once by assigning the function.
    # Both functions share the same arguments.
    if config.reuse_cache:
        primary_source = cached_distance_first
    else:
        primary_source = get_table_distance

    for s in subject_mappings:
        result = primary_source(config, query_mapping.species_id, s.species_id)
        if default is not None and result is None:
            yield default
        else:
            yield result


class species_pair_distance:
    """ A small class to store the IDs of a species pair and their ML protein
    distance. """
    __slots__ = ['species1', 'species2', 'distance']

    def __init__(self, spec1, spec2, dist):
        self.species1 = spec1
        self.species2 = spec2
        self.distance = dist

    @property
    def spec_set(self):
        """ Get both species as a species set to compare identity regardless
        of order. """
        # We do not need or want this set to be mutable.
        return frozenset([self.species1, self.species2])

    def equals(self, other):
        return self.spec_set == other.spec_set


if __name__ == "__main__":
    """ This script is only executed as a script for testing purposes! """
    config_file = Path('test.config')
    config = set_params(config_file)

    print('Testing the integrity of the species mapping file!')

    test_mapping = species_mapping(config)
    print(test_mapping)
