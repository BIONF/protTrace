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

import os
import sys
import time
import argparse
from pathlib import Path

from multiprocessing.pool import Pool
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.Align import MultipleSeqAlignment as msa
from Bio.Phylo.Consensus import bootstrap

from utils.configure import set_params
from utils.data_api import gen_proteomes, gen_pairwise_orthologs

# Supporting script for maximum likelihood calculations #


def calculate_species_distances(config):
    """ Orchestrates the calculation of species distances between
    the query species and all other species in the
    speciesTreeMapping list. <config> is expected as an
    instance of the ProtTrace configuration (generated with
    configure.setParams) """

    # Read the query species from the ProtTrace configuration
    # The special species, 'ALL', updates all species in the mapping file
    query = config.species
    species_mapping = config.species_map

    def read_species_set(species_mapping_file):
        """ Gather the set of all subject species in this project
        from the speciesTreeMapping file. """
        species_set = set()
        try:
            with open(species_mapping_file, 'r') as sm:
                for line in sm:
                    species_set.add(line.split()[-1])
            return species_set
        except IOError:
            sys.exit(print_error('Could not open {0}. Please check the path.'
                     .format(species_mapping_file)))

    prepare_directories(config)
    previous_work_dir = move_working_dir(None, config)

    def search_and_complement_missing_species(query, species_set, config,
                                              check_table=True,
                                              check_cache=True):

        # Set the default set of missing species
        missing_species = species_set

        def check_table_for_missing_species(missing_species,
                                            query, species_set, config):

            if check_table and len(missing_species) > 0:
                # Exclude all species with existing distances in
                # species_max_likelihood. The function will return a SET of
                # species. Their order when iterated will be random!
                missing_species = check_species_max_likelihood_table(
                    query,
                    species_set,
                    config)

            return missing_species

        def check_cache_for_missing_species(missing_species, query, config):

            if check_cache and len(missing_species) > 0:
                if not os.path.exists(config.path_cache):
                    return missing_species
                # Exclude missing species with existing .lik files in the cache
                # directory. We check the cache directory for
                # {query species}_{missing species}.lik files
                missing_species = {species for species in missing_species
                                   if not os.path.isfile("{0}/{1}_{2}.lik"
                                                         .format(
                                                            config.path_cache,
                                                            query, species))}
            return missing_species

        missing_species = check_table_for_missing_species(missing_species,
                                                          query,
                                                          species_set,
                                                          config)
        missing_species = check_cache_for_missing_species(missing_species,
                                                          query,
                                                          config)

        # If we are still missing species distances, we calculate them now
        # All distances are copied to the cache directory. This means they stay
        # backed up when the not saving cache option would wipe them
        if len(missing_species) > 0:

            if config.nr_processors > 1 and len(missing_species) > 1:
                # Every pair of species is parallelized. Each pairwise distance
                # computation gets one core.
                with Pool(config.nr_processors) as species_pair_process_pool:
                    species_pair_process_pool.map(
                        calculate_protein_distances_parallelized,
                        [[query, species, config, config.path_cache, 1, 1, 0]
                         for species in missing_species
                         if species is not query])

            # One species pair is processed at once. This means, the function
            # is allowed to use the full number of available cores
            # (by passing None)
            else:
                for species in missing_species:
                    if species is not query:
                        calculate_protein_distances(query, species, config,
                                                    config.path_cache,
                                                    None, 0)
    # DEBUG
    # The species pair YEAST-SULMD only has around 86 pairwise orthologs in OMA
    # database, release August 2020.
    # search_and_complement_missing_species("YEAST", set(["SULMD"]), config,
    # False, False)
    # sys.exit()

    species_set = read_species_set(species_mapping)

    if query == 'ALL':
        # Here, we use each species as a query
        for species in species_set:
            # Existing distances in the cache are excluded, but this routine
            # replaces distances in the ML table
            search_and_complement_missing_species(species, species_set, config,
                                                  False, True)
    else:
        search_and_complement_missing_species(query, species_set, config)

    # Return to the previous working directory
    move_working_dir(previous_work_dir, config)


def prepare_directories(config, target_dir=None, species1=None, species2=None):
    """ Creates missing directories which are necessary later. """

    if not os.path.exists(config.path_distance_work_dir):
        os.mkdir(config.path_distance_work_dir)

    if target_dir is not None:
        if not os.path.exists(config.path_cache):
            os.mkdir(config.path_cache)

    if species1 is not None and species2 is not None:
        species_work_dir = '{0}/{1}_{2}'.format(
            config.path_distance_work_dir, species1, species2)
        if not os.path.exists(species_work_dir):
            os.mkdir(species_work_dir)
        return species_work_dir


def move_working_dir(previous_work_dir, config, species_dir=None):
    """ Dynamically moves the current directory based on the intended
    workflow. """

#    def unit_test_control(expected_previous, expected_current):
#        print(previous_work_dir)
#        print(os.getcwd())
#
#    # This will be properly implemented with assert when the unit test update
#    # comes. Now, it just prints the current state
#    unit_test_control(None, None)

    # Species_dir is specified
    if species_dir is not None:
        # species_dir -> distance_dir
        if os.getcwd() == species_dir:
            os.chdir(previous_work_dir)
            return species_dir
        else:
            # distance_dir -> species_dir
            if os.getcwd() != config.path_distance_work_dir:
                os.chdir(config.path_distance_work_dir)
            os.chdir(species_dir)
            return config.path_distance_work_dir

    # distance_dir -> prottrace_work_dir
    # distance_dir -> anywhere if executed as __main__
    if os.getcwd() == config.path_distance_work_dir:
        os.chdir(previous_work_dir)
        return config.path_distance_work_dir

    # Move to the dedicated distance calculation root directory for all species
    # after noticing the previous directory to move back later
    # prottrace_work_dir -> distance_dir
    # anywhere -> distance_dir if executed as __main__
    previous_work_dir = os.getcwd()
    if previous_work_dir != config.path_distance_work_dir:
        os.chdir(config.path_distance_work_dir)
        return previous_work_dir

    print_error('The current working directory could not be recognized. '
                'Exiting.')
    sys.exit()


def delete_temporary_files(species_1, species_2, fa_dir=None, aln_dir=None,
                           bootstrap_count=0, result_file=None):

    # If there is a result file, but not completed, restrain from
    # deleting temporary files for reuse
    if result_file is not None:
        if (not os.path.exists(result_file)
                or len(open(result_file).read().split('\n')) == 0):
            return
    try:
        os.system('rm -rf temp_puzzleParams.txt')
        os.system(f'rm -rf {species_1}_{species_2}.phy*')
        # This check rather exists to prevent file not found errors
        if bootstrap_count > 0:
            os.system('rm -rf {species_1}_{species_2}_bootstrap_*')
        print_progress('Cleanup finished.')
    except FileNotFoundError:
        pass


# Checks the species maximum likelihood file for existing distances
# Returns missing target species OMA IDs
def check_species_max_likelihood_table(query, targets, config):
    """ Checks the species maximum likelihood file for existing distances
    and returns species IDs with missing valid distace numbers. """

    # Check whether all query-subject pairs have already
    # computed distances in the species_max_likelihood file
    species_max_likelihood = config.species_MaxLikMatrix
    # We need a list first, because we first fetch the names
    # in the first line, then we fetch existing numbers
    # in the query species row, where both indices must
    # coalign
    computed_species_list = list()
    computed_species_set = set()

    def generate_max_likelihood_distance_table_rows(filereader):
        for row in filereader:
            yield row.rstrip().split('\t')

    try:
        with open(species_max_likelihood, 'r') as ml:
            # The first row contains all species names
            # Its first column is empty
            # We populate the species list with it
            # The file reader will be positioned at the next line anyways
            computed_species_list = ml.readline().rstrip().split("\t")[1:]
            # Now we look up our query species row for which distances are
            # actually computed / not N/A / not missing
            for split_row in generate_max_likelihood_distance_table_rows(ml):
                # Find the row where the first item corresponds to our query
                if split_row[0] == query:
                    # The first item was only necessary to recognize the query
                    # This way, we save a +1 index operation for every column
                    # (and the headache)
                    split_row = split_row[1:]
                    # Create a set of all species names whose index in this
                    # column contains a valid digit
                    computed_species_set = {
                        computed_species_list[i] for i in range(len(split_row))
                        if split_row[i] is not None and split_row[i].isdigit()}
                    break
    except IOError:
        sys.exit(print_error('Could not open {0}. Please check the path.'
                 .format(species_max_likelihood)))

    # This is done to remove references and in case
    # targets is a list
    target_species_set = set(targets)

    # Returns species in the speciesTreeMapping list without
    # digits in SpeciesMaxLikelihood
    return target_species_set.difference(computed_species_set)


def gen_bootstraps(args):
    """ Converts the concatenated alignment SeqRecords into a MSA for
    Bio.Phylo.Consensus.bootstrap to work with. Then returns the
    bootstrap SeqRecords. """

    try:
        align_records = args[0]
        bootstrap_count = args[1]

        bootstraps = bootstrap(msa(align_records),
                               bootstrap_count)

        # We unpack the list of bootstraps. Then we unpack the
        # MultipleSeqAlignment objects into a List of SeqRecord objects by
        # calling the MultipleSeqAlignment iterator. A list comprehension
        # somehow always yields a MultipleSeqAlignment object, no matter how
        # many times of iteration.
        records = []
        for bs in bootstraps:
            for record_list in bs:
                records.append(record_list)

        return records
    except Exception as e:
        raise e


def concatenate_alignment(species_1, species_2, alignment_generator,
                          bootstrap_count=0, nr_processors=1):

    concat_1 = ''
    concat_2 = ''
    starttime = time.time()
    count = 0
    for alignment in alignment_generator:
        concat_1 += alignment[0]
        concat_2 += alignment[1]
        # The progress is counted when a new batch of pairwise protein
        # alignments are called.
        if count % 1000 == 0:
            print_progress(f'Aligned {count} sequence pairs '
                           'over {:.0f} seconds'.format(time.time()
                                                        - starttime),
                           guarantee_print=True)
        count += 1

    concat_alignments = (SeqRecord(Seq(concat_1), id=species_1),
                         SeqRecord(Seq(concat_2), id=species_2))

    def duplicate_seqs(species_1, species_2, records):
        yield records[0]
        yield SeqRecord(records[0].seq, id=records[0].id + '_dub')
        yield records[1]
        yield SeqRecord(records[1].seq, id=records[1].id + '_dub')

    # Perform a bootstrap analysis on the alignment and note the created
    # variance
    # def gen_bootstraps_legacy(alignment_length, generated_bootstrap_count):

    #     # Seed the singleton random module using the current system clock
    #     # time (by passing no parameter)
    #     seed = random.randrange(sys.maxsize)
    #     random.seed(seed)

    #     print('The seed for bootstrapping the alignment is: ', seed)
    #     with open('bootstrap_sample_rng_seed.R', 'w') as seed_file:
    #         seed_file.write('# The seed for sampling columns for pairwise '
    #                         'distance\n# bootstrap analysis is:\n# '
    #                         + str(seed) + '\n')

    #     # Generate a list of sequence indices to sample
    #     # The list is sampled from all alignment indices with equal weights
    #     # and replacement.
    #     for c in range(generated_bootstrap_count):
    #         yield random.choices(range(alignment_length), k=alignment_length)

    def write_phylip(spec_1, spec_2, concats, suffix=''):
        SeqIO.write(duplicate_seqs(spec_1, spec_2, concats),
                    f'{species_1}_{species_2}{suffix}.phy', 'phylip')

    # Write the main concatenated alignment into a phylip file to calculate
    # the species distance.
    write_phylip(species_1, species_2, concat_alignments)

    if bootstrap_count > 0:
        count = 0
        args = [(concat_alignments, 1) for b in range(bootstrap_count)]
        with Pool(nr_processors) as process_pool:
            for sample in process_pool.imap_unordered(gen_bootstraps, args):
                count += 1
                write_phylip(species_1, species_2, sample,
                             suffix=f'_bootstrap_{count}')


def search_protein_pairs(species_1, species_2, config):
    """ Searches the given table for mutual presences of
    species 1 and species 2 within the first two columns
    and returns both columns verbatim in the order the
    species were passed to this function. """

    def generate_large_orthologous_pair_file_buffers(pairs_src):
        # Read in oma_pair lines that contain proteins
        # of both species together
        with open(pairs_src, 'r') as pair_mapping:
            tmp_lines = pair_mapping.readlines(134217728)
            while tmp_lines:
                # Identify the input species pair
                # among the orthologous protein pairs
                yield tmp_lines
                startingtime = time.time()
                tmp_lines = pair_mapping.readlines(134217728)
                print('Read: {0}'.format(str(time.time() - startingtime)))

    def gen_pairwise_orth_lines(species_1, species_2, pairs_src):
        with open(pairs_src, 'r') as pair_table:
            for line in pair_table:
                if species_1 in line and species_2 in line:
                    yield line

    # def read_pairs_src(pairs_src, src_type):
    #     if src_type == 'oma':
    #         src_func = pairs_src.get_pairwise_orthologs
    #         src_args = [species_1, species_2]
    #     elif src_type == 'file':
    #         src_func = gen_pairwise_orth_lines
    #         src_args = [species_1, species_2, pairs_src]

    #     for line in src_func(*src_args):
    #         yield line

    def order_columns(columns, species_1, species_2, src_type):
        if src_type == 'oma':
            col_1 = 2
            data_col_1 = 0
            col_2 = 3
            data_col_2 = 1
        elif src_type == 'file':
            col_1, data_col_1 = 0
            col_2, data_col_2 = 1

        if species_1 in columns[col_1]:
            return (columns[data_col_1], columns[data_col_2])
        else:
            return (columns[data_col_2], columns[data_col_1])

    try:
        sequence_count = 1
        src_type = 'oma'
        for line in gen_pairwise_orthologs(config, species_1, species_2):
            columns = line.rstrip().split('\t')
            # The previous generator already established that
            # both species are present in the columns. Here, we
            # look at the exact column position.
            yield order_columns(columns, species_1, species_2, src_type)
            sequence_count += 1
            # Inform the user about the current progress
            if sequence_count % 1000 == 0:
                print_progress('Sequence Nr.: {0}'.format(str(sequence_count)))
        print_progress('Sequence Nr.: {0}'.format(str(sequence_count)))

        if sequence_count == 1:
            print_warning('No pairwise orthologs found!')
            return
    except KeyboardInterrupt:
        sys.exit('The user interrupted the compilation of orthologous pairs!')
    except FileNotFoundError:
        sys.exit(print_error('The OMA pairs file is missing!'))


def provide_protein_pair_sequences(species1, species2, sequence_source):
    """ Provide a dict of protein ids and sequences. """

    try:
        return sequence_source.index_multiple_seqs([species1, species2])
        # return sequence_source.index_seqs()
    except KeyboardInterrupt:
        print('Gathering the sequences of pairwise orthologs was interrupted '
              'by the user.')
        sys.exit()
    except FileNotFoundError:
        print_error('The FASTA file that contains the sequences is missing!')


def align_pairwise2_parallelized(args):
    """ Provides a top-level function to parallelize the alignment of two
        sequences, args[0] and args[1]. Other parameters are static for every
        pair of sequences."""

    # one_alignment_only retrieves the first of the equally best scoring
    # alignments.
    # d: match scores are derived from a list (the substitution matrix)
    # s: same gap penalties for both sequences (values for gapopen and
    # gapextend were taken from MAFFT 6 L-INS-i
    # One alignment consists of a tuple of four elements:
    # sequence 1, sequence 2, score, beginning, 0-based character count

    starttime = time.time()
    alignment = pairwise2.align.globalds(args[0], args[1],
                                         substitution_matrices.load(
                                             'BLOSUM62'), -10.0, -1.0,
                                         penalize_end_gaps=True,
                                         one_alignment_only=True)[0]

    print('Aligned: {0}'.format(str(time.time() - starttime)), flush=True)

    # Somewhere between Bio 1.7.2 and 1.7.8, pairwise2 now returns an
    # Alignment object with seqA and seqB, instead of a list containing
    # the sequences in element 0 and element 1. The module object "Alignment"
    # cannot be pickled back.
    return (str(alignment.seqA), str(alignment.seqB))


def align_PairwiseAlignment_parallelized(args):
    """ Provides a top-level function to parallelize calls to PairwiseAlignment
    .aligner (args[0]). The two sequences are args[1] and args[2]. """

    try:
        # The PairwiseAligner object stores the settings of the alignment
        # process. This object cannot be pickled, so it has to be created in
        # every task again.
        aligner = Align.PairwiseAligner(open_gap_score=-10.0,
                                        extend_gap_score=-1.0)
        aligner.mode = 'global'
        aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')

        # starttime = time.time()
        # We remove the alignment indicator bars by removing the second line.
        alignment = str(aligner.align(args[0], args[1])[0]).split('\n')[0:3:2]
        # print('Aligned: {0}'.format(str(time.time() - starttime)))
        return (alignment[0], alignment[1])
    except Exception as e:
        print_error('A problem has occurred when aligning the proteins!')
        print(args)
        raise e


def calculate_protein_distances_parallelized(args):
    """ A pickle-able dispatcher for calculate_protein_distances. """
    calculate_protein_distances(*args)


# Calculates the pairwise species maximum likelihood distance
# between species1 and species2. The distance is copied to the
# target_dir. The config is loaded from ProtTrace's configure.py
# if this script is executed directly.
def calculate_protein_distances(species1, species2,
                                config, target_dir,
                                preset_nr_processors=None,
                                bootstrap_count=0):

    print_progress('Calculating the protein distance between {0} and {1}.'
                   .format(species1, species2))

    treepuzzle = config.treepuzzle
    if treepuzzle is None:
        print_error('No path to TreePUZZLE has been configured!')
        sys.exit()
    if preset_nr_processors is None:
        nr_processors = config.nr_processors
    else:
        nr_processors = preset_nr_processors
    delete_temp = config.delete_temp

    """ Directory management """

    work_dir = prepare_directories(config, target_dir, species1, species2)
    root_dir = move_working_dir(None, config, work_dir)

    print_progress('Aligning and concatenating pairwise orthologous protein '
                   'sequences.')

    def generate_PairwiseAlignment_arguments(prot_pairs_generator,
                                             proteomes):
        """ Collects the sequences of pairwise orthologous proteins. """
        for prot_pair in prot_pairs_generator:
            yield (str(proteomes[0].get_protein(prot_pair[0]).seq),
                   str(proteomes[1].get_protein(prot_pair[1]).seq))

    def align_pairwise_proteins(prot_pairs_generator,
                                proteomes,
                                nr_processors):
        """ Aligns all pairwise orthologous protein sequences between species1
            and species2 using the Biopython pairwise2 module. If the distance
            calculation between species is done sequentially, then this process
            can spawn child processes by providing a nr_processors count
            greater than 1. """

        if nr_processors > 1:
            try:
                with Pool(nr_processors) as align_processes:
                    for alignment in align_processes.imap_unordered(
                            align_PairwiseAlignment_parallelized,
                            generate_PairwiseAlignment_arguments(
                            prot_pairs_generator, proteomes),
                            1000):
                        yield alignment

            except KeyboardInterrupt:
                sys.exit('The user interrupted the alignment of pairwise '
                         'orthologous proteins.')

        else:
            try:
                for arguments in generate_PairwiseAlignment_arguments(
                        prot_pairs_generator, proteomes):
                    yield align_PairwiseAlignment_parallelized(arguments)

            except KeyboardInterrupt:
                sys.exit('The user interrupted the alignment of pairwise'
                         'orthologous proteins.')
            # Since this code section will be executed if the species pairs are
            # parallelized, we need to explicity propagate any exception ot the
            # parent process.
            except Exception as e:
                os.chdir(root_dir)
                raise e

    """ Load proteome sequences. """

    proteomes = list(gen_proteomes(config, species1, species2))

    """ Align pairwise orthologous sequences and concatenate them. """

    # The collection of orthologous pairs, their sequences, their alignments
    # and their concatenation are all chained with generators here.
    concatenate_alignment(species1, species2,
                          align_pairwise_proteins(
                              search_protein_pairs(species1, species2, config),
                              proteomes, nr_processors),
                          bootstrap_count, nr_processors)

    """ Calculate the pairwise species distance. """

    # Executes TreePUZZLE to calculate the distance between the species pair
    def calculate_pairwise_distance(concat_file, treepuzzle):

        # Prepare the parameter file for TreePUZZLE
        with Path('temp_puzzleParams.txt').open('w') as pp:
            pp.write(str(concat_file) + '\ne\nm\nm\nm\nm\nm\nm\ny\n')

        # Execute TreePUZZLE
        os.system(f'{treepuzzle} < temp_puzzleParams.txt >/dev/null')

    def read_calculated_distance(concat_file):
        try:
            with concat_file.with_suffix('.phy.dist').open('r') as concat:
                next(concat)
                return concat.readline().split()[3]
        except FileNotFoundError:
            print_error('Error: The pairwise distance calculation result file '
                        'from TreePUZZLE cannot be found!')
            sys.exit()

    def write_result(distance, filename):
        with filename.open('w') as output:
            output.write(distance + '\n')

    def write_main_result(distance, filename, target_dir):

        write_result(distance, filename)

        # Copy the main output file to the given target directory, if given
        # The target directory is originally designed to be the ProtTrace cache
        # directory
        if target_dir is not None and target_dir.is_dir():
            write_result(distance, target_dir / filename)

    # Process the pairwise species maximum likelihood distance
    print_progress('Calculating and writing the maximum likelihood distance!')
    concat_filename = Path(f'{species1}_{species2}.phy')
    result_file = concat_filename.with_suffix('.lik')
    calculate_pairwise_distance(concat_filename, treepuzzle)
    main_distance = read_calculated_distance(concat_filename)
    write_main_result(main_distance, result_file, target_dir)

    # Process the distances calculated from bootstrapped alignments
    # First, we define a couple of generator functions

    # Anticipate the output file names of the
    # bootstrapped concatenated alignments. For aesthetic reason, they start at
    # 1.
    def generate_bootstrap_concat_filenames(bootstrap_count):
        for i in range(1, bootstrap_count + 1):
            yield Path(f'{species1}_{species2}_bootstrap_{i}.phy')

    def generate_bootstrap_distances(bootstrap_count, treepuzzle):

        for bootstrap_concat_filename in generate_bootstrap_concat_filenames(
                bootstrap_count):
            calculate_pairwise_distance(bootstrap_concat_filename, treepuzzle)
            yield read_calculated_distance(bootstrap_concat_filename)

    def calculate_and_write_bootstrap_distance_table_output(main_distance,
                                                            bootstrap_count,
                                                            treepuzzle):

        # Write the main distance and its bootstrap distances into a table
        # for followup statistical analyses
        # The bootstrap column is a R friendly boolean to separate the main
        # distance from the bootstrap values
        separator = '\t'
        columns = [separator.join(['Species_1', 'Species_2',
                                   'Distance', 'Bootstrap'])]
        columns.append(separator.join([species1, species2,
                                       main_distance, 'FALSE']))
        for distance in generate_bootstrap_distances(bootstrap_count,
                                                     treepuzzle):
            columns.append(separator.join([species1, species2,
                                           distance, 'TRUE']))
        with open(species1 + '_' + species2 + '_bootstrap.lik', 'w')\
                as analysis_output:
            analysis_output.write('\n'.join(columns) + '\n')

    if bootstrap_count > 0:
        calculate_and_write_bootstrap_distance_table_output(main_distance,
                                                            bootstrap_count,
                                                            treepuzzle)

    if delete_temp:
        delete_temporary_files(species1, species2,
                               bootstrap_count, result_file)

    print_progress(f'Finished species pair {species1} - {species2}',
                   guarantee_print=True)
    move_working_dir(root_dir, config)


def print_progress(message, guarantee_print=False):
    print('### ' + message + ' ###', flush=guarantee_print)


def print_warning(message):
    print('Warning: ' + message)


def print_error(message):
    print('ERROR: ' + message)


def main():
    """ The entry point to calculate the distance between two species. """

    def Argparse():
        """ Parses the arguments into an object. """

        parser = argparse.ArgumentParser(description='Calculates the maximum '
                                         'lokelihood protein distance between '
                                         'two species')
        parser.add_argument('-q', '--query', type=str, nargs='?', default=None,
                            help='The query species. It is taken from the '
                            'config per default.')
        parser.add_argument('-t', '--target', type=str, help='The target '
                            'species.')
        parser.add_argument('-c', '--config', type=str, help='The path to the '
                            'configuration file.')
        parser.add_argument('-p', '--processors', type=int, default=None,
                            help='Number of cores assigned.')
        parser.add_argument('-b', '--bootstraps', type=int, default=0,
                            help='Number of bootstrap alignments')
        return parser.parse_args()

    arguments = Argparse()

    query = arguments.query
    target = arguments.target
    config_file = os.path.abspath(arguments.config)
    nr_processors = arguments.processors

    config = set_params(config_file)

    if query is None:
        query = config.species

    target_dir = None

    if os.path.exists(config.path_cache):
        target_dir = Path(config.path_cache)

        # Check whether any combination of both species exist.
        # Their distances should be identical.
        precomputed_1 = Path(f'{target_dir}/{query}_{target}.lik')
        precomputed_2 = Path(f'{target_dir}/{target}_{query}.lik')
        if precomputed_1.exists() or precomputed_2.exists():
            print('Species distance is already computed!')
            return 0

    # If this script is started on its own calculate the pairwise distance
    # between the specified species pair
    calculate_protein_distances(arguments.query, arguments.target, config,
                                target_dir, preset_nr_processors=nr_processors,
                                bootstrap_count=arguments.bootstraps)

    return 0


# This defines the start of this script if someone wants to calculate
# distances separately with this script
if __name__ == "__main__":
    sys.exit(main())
