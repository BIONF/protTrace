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

import os, sys
import glob
import random
import time
import subprocess
from multiprocessing.pool import ThreadPool
from multiprocessing.pool import Pool
from Bio import SeqIO
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist

### Supporting script for maximum likelihood calculations ###

def calculate_species_distances(config):
    """ Orchestrates the calculation of species distances between
    the query species and all other species in the
    speciesTreeMapping list. <config> is expected as an
    instance of the ProtTrace configuration (generated with
    configure.setParams) """

    # Read the query species from the ProtTrace configuration
    # The special species, 'ALL', updates all species in the mapping file
    query = config.species
    species_mapping = config.hamstr_oma_tree_map

    # Gather the set of all subject species in this project
    # from the speciesTreeMapping file
    species_set = set()
    try:
        with open(species_mapping, 'r') as sm:
            for line in sm:
                species_set.add(line.split()[-1])
    except IOError:
        sys.exit('ERROR: Could not open %s. Please check the path.' %crossRefFile)

    # Move to the dedicated distance calculation root
    # directory for all species after noticing the previous directory to move back later
    previous_work_dir = os.getcwd()
    os.chdir(config.path_distance_work_dir)

    # Check and prepare the cache location
    if not os.path.exists(config.path_cache):
        os.mkdir(config.path_cache)

    def search_and_complement_missing_species(query, species_set, config, check_table=True, check_cache=True):

        # Set the default set of missing species
        missing_species = species_set

        if check_table and len(missing_species) > 0:
            # Exclude all species with existing distances in species_max_likelihood
            # The function will return a SET of species. Their order when iterated will be random!
            missing_species = check_species_max_likelihood_table(query,species_set,config)

        if check_cache and len(missing_species) > 0:
            # Exclude missing species with existing .lik files in the cache directory
            # We check the cache directory for {query species}_{missing species}.lik files
            missing_species = {species for species in missing_species if not os.path.isfile("{0}/{1}_{2}.lik".format(config.path_cache,query,species))}

        # If we are still missing species distances, we calculate them now
        # All distances are copied to the cache directory
        # This means they stay backed up when the not saving cache option would wipe them
        if len(missing_species) > 0:

            # The first, commented parallelization solution is only possible if hierarchical
            # multiprocessing is implemented, i.e. child processes are allowed to spawn child
            # processes on their own. Unless that is ensured, the second parallelized solution
            # can be used, where each species pair is computed in parallel.

            # These two questions ensure that the load is worth the parallelization
            #if config.nr_processors >= 4 and len(missing_species) > (config.nr_processors // 4):
            #if True:
            #    # If more than 7 cores are available, the number of cores is split by 4 to
            #    # parallelize processing
            #    with Pool(config.nr_processors // 4) as species_pair_process_pool:
            #        species_pair_process_pool.map(calculate_protein_distances_parallelized, [[query,species,config,config.path_cache,4,1,0] for species in missing_species if species is not query])
            if config.nr_processors > 1 and len(missing_species) > 1:
                # Every pair of species is parallelized. Each pairwise distance computation
                # gets one core.
                with Pool(config.nr_processors) as species_pair_process_pool:
                    species_pair_process_pool.map(calculate_protein_distances_parallelized, [[query,species,config,config.path_cache,1,1,0] for species in missing_species if species is not query])

            # One species pair is processed at once. This means, the function is allowed
            # to use the full number of available cores (by passing None)
            else:
                for species in missing_species:
                    if species is not query:
                        calculate_protein_distances(query,species,config,config.path_cache,None,1,0)
    #DEBUG
    #search_and_complement_missing_species("YEAST",set(["SULMD"]),config,False,False)
    #sys.exit()

    if query == 'ALL':
        # Here, we use each species as a query
        for species in species_set:
            # Existing distances in the cache are excluded, but this routine
            # replaces distances in the ML table
            search_and_complement_missing_species(species,species_set,config,False,True)
    else:
        search_and_complement_missing_species(query,species_set,config)

    # Return to the previous working directory
    os.chdir(previous_work_dir)

def delete_temporary_files(species1, species2, fa_dir, aln_dir, bootstrap_count = 0, result_file = None):
    # If there is a result file, but not completed, restrain from
    # deleting temporary files for reuse
    if result_file is not None:
        if not os.path.exists(result_file) or len(open(result_file).read().split('\n')) == 0:
            return
    try:
        os.system('rm -rf {0} {1} temp_puzzleParams.txt'.format(fa_dir, aln_dir))
        os.system('rm -rf {0}_{1}.phy*'.format(species1, species2))
        # This check rather exists to prevent file not found errors
        if bootstrap_count > 0:
            os.system('rm -rf {0}_{1}_bootstrap_*'.format(species1, species2))
        print('Cleanup finished.')
    except FileNotFoundError:
        pass

# Checks the species maximum likelihood file for existing distances
# Returns missing target species OMA IDs
def check_species_max_likelihood_table(query,targets,config):
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
                # Find the row where the first item corresponds to our query species
                if split_row[0] == query:
                    # The first item was only necessary to recognize the query species
                    # This way, we save a +1 index operation for every column
                    # (and the headache)
                    split_row = split_row[1:]
                    # Create a set of all species names whose index in this column
                    # contains a valid digit
                    computed_species_set = {computed_species_list[i] for i in range(len(split_row)) if split_row[i] is not None and split_row[i].isdigit()}
                    break
    except IOError:
        sys.exit('ERROR: Could not open %s. Please check the path.' %crossRefFile)

    # This is done to remove references and in case
    # targets is a list
    target_species_set = set(targets)

    # Returns species in the speciesTreeMapping list without
    # digits in SpeciesMaxLikelihood
    return target_species_set.difference(computed_species_set)

# The multiprocessed function needs to be pickle-abled
# Therefore, it must be defined at top-level
def sequence_pair_to_phylip_multiprocessed(argument_list):
    """ A pickle-able dispatcher for sequence_pair_to_phylip. """
    sequence_pair_to_phylip(*argument_list)

# This function is multiprocessed for bootstrapped alignments
# Therefore, it must be defined at top-level
#def sequence_pair_to_phylip(species1,species2,sequences,subalignment_positions,sampled_indices=None,bootstrap_index=None):
def sequence_pair_to_phylip(species1,species2,sequences,sampled_indices=None,bootstrap_index=None):

    # Here, we compile the bootstrapped sequences to continuous strings
    # The regular concatenated protein alignments are just made continuous
    # Both species remain separated
    local_sequences = ["",""]
    if sampled_indices is not None:
        if len(sampled_indices) != len(sequences[0]):
            print("ERROR: The bootstrap sequence has not the same length as the protein pair count!")
        # Sample the amino acids from both aligned sequences by their
        # common index in the alignment
        for sampled_index in sampled_indices:
            local_sequences[0] += sequences[0][sampled_index]
            local_sequences[1] += sequences[1][sampled_index]
    else:
        local_sequences = sequences

    def append_phylip_seq_with_spaces(begin,seq):
        # The sequence is indented to the right by 10 spaces
        appended_seq = " " * 10
        # The next 50 positions are inserted with one space between
        # 10 positions. The minimum function ensures that we do not
        # try to access empty sequence positions
        for i in range(begin,begin + min(len(seq) - begin,50),10):
            appended_seq += " " + seq[i:i+10]
            # Fill an open ending block with hyphens
            appended_seq += "-" * (10 - len(seq[i:i+10]))
        appended_seq += "\n"
        return appended_seq

    def append_phylip_seq_with_spaces_with_species(species,seq):
        appended_seq = species + " " * (10 - len(species))
        for i in range(0,50,10):
            appended_seq += " " + seq[i:i+10]
        appended_seq += "\n"
        return appended_seq
    # This function produces a PHYLIP file where the full sequence is written
    # into one line. By commenting out the for loop for filling the interleaved
    # format, the alignment can be bootstrapped with SEQBOOT for confirmation.
    #def append_phylip_seq_with_spaces_with_species_singleline(species,seq):
    #    appended_seq = species + " " * (10 - len(species))
    #    for i in range(0,len(seq),10):
    #        appended_seq += " " + seq[i:i+10]
    #    appended_seq += "\n"
    #    return appended_seq

    # Build the PHYLIP file from the continuous alignment
    alignment = ""
    # The first line contains the entire length of the original protein sequence
    # Here, we just copy the first and last position
    alignment += "4 " + str(len(local_sequences[0])) + "\n"
    # The next 2 lines contain the respective species ids of each line
    # We duplicate the sequence to generate an imaginative
    # 4 species alignment. We need this to calculate a
    # tree and calculate the distance between sister
    # clades (each clade is the same species twice)
    alignment += append_phylip_seq_with_spaces_with_species(species1,local_sequences[0])
    alignment += append_phylip_seq_with_spaces_with_species(species1 + "_dub",local_sequences[0])
    alignment += append_phylip_seq_with_spaces_with_species(species2,local_sequences[1])
    alignment += append_phylip_seq_with_spaces_with_species(species2 + "_dub",local_sequences[1])

    # Adds spacing between alignment blocks
    alignment += "\n"
    # This will fill the rest of the alignment file. We assume that both sequences
    # are equally long
    for i in range(50, len(local_sequences[0]), 50):
        alignment += append_phylip_seq_with_spaces(i,local_sequences[0])
        alignment += append_phylip_seq_with_spaces(i,local_sequences[0])
        alignment += append_phylip_seq_with_spaces(i,local_sequences[1])
        alignment += append_phylip_seq_with_spaces(i,local_sequences[1])
        # Adds spacing between alignment blocks
        alignment += "\n"

    # Assemble the name of the concatenation file
    # The bootstrap concatenation file is named as such
    concat_file_name = species1 + '_' + species2
    if sampled_indices is not None:
        concat_file_name += '_bootstrap_{0}'.format(bootstrap_index)
    concat_file_name += '.phy'

    with open(concat_file_name, 'w') as concat_file:
        concat_file.write(alignment)
    # Only the non-bootstrap result file receives a
    # subalignment_positions file
    #if sampled_indices is None:
    #    with open('subalignment_positions.txt', 'w') as subpositions_file:
    #        for line in subalignment_positions:
    #            # Columns are separated with spaces
    #            subpositions_file.write(' '.join([str(l) for l in line]) + "\n")

def concatenate_alignment(species1, species2, alignment_generator, bootstrap_count = 0, nr_processors = 1):

    concatenated_alignments = ["",""]
    for alignment in alignment_generator:
        concatenated_alignments[0] += alignment[0]
        concatenated_alignments[1] += alignment[1]

    # Perform a bootstrap analysis on the alignment and note the created variance
    def generate_bootstraps(alignment_length, generated_bootstrap_count):

        # Seed the singleton random module using the current system clock time
        # (by passing no parameter)
        seed = random.randrange(sys.maxsize)
        random.seed(seed)

        print('The seed for bootstrapping the alignment is: ', seed)
        with open('bootstrap_sample_rng_seed.R','w') as seed_file:
            seed_file.write('# The seed for sampling columns for pairwise distance\n# bootstrap analysis is a follows:\n# '+ str(seed) + '\n')

        # Generate a list of sequence indices to sample
        # The list is sampled from all alignment indices with equal weights and replacement
        for c in range(generated_bootstrap_count):
            yield random.choices(range(alignment_length), k=alignment_length)

#    if bootstrap_count > 0:
#        bootstrap_alignment_indices = create_bootstraps(len(sequences[0]),bootstrap_count)

    print('Concatenate and duplicate the pairwise aligned sequences!')
    # Generate the main concatenated and duplicated alignment for calculating
    # the pairwise species distance
    sequence_pair_to_phylip(species1,species2,concatenated_alignments)

    def generate_multiprocessing_to_phylip_args_list(species1, species2, sequences, subalignment_positions, bootstrap_count):
        i = -1
        # Sample a set of position indices to draw aligned amino acids from both sequences
        for bootstrap_position_indices in generate_bootstraps(len(sequences[0]), bootstrap_count):
            i += 1
            yield [species1,species2,sequences,subalignment_positions,bootstrap_position_indices,i]

    if bootstrap_count > 0:
        print('Concatenate and duplicate the bootstrapped pairwise alignments!')
        # Generate the bootstrapped versions of the concatenated alignment for
        # estimating the distance variance
        # Each alignment bootstrap creates a separate process
        # The pool iterates through "bootstrap_alignment_indices"
        try:
            with Pool(nr_processors) as process_pool:
                # imap is a lazy loader of the bootstrap position indices
                # Sampling every position and putting their integers in a list takes
                # much more RAM than when the positions are written into two strings
                process_pool.imap(sequence_pair_to_phylip_multiprocessed, generate_multiprocessing_to_phylip_args_list(species1,species2,sequences,subalignment_positions,bootstrap_count),5)
        except KeyboardInterrupt:
            sys.exit('The user interrupted the generation of bootstrapped alignments!')

def concatenate_alignment_legacy(species1, species2, alignment_count, alignment_directory, bootstrap_count = 0, nr_processors = 1):

    # The sequences of both species are separated for better handling
    # For bootstrapping, the first index tells us the aligned protein pair
    # The second index tells us the species
    sequences = ["",""]
    protein_pair = 0
    current_species = 2
    subalignment_positions = []
    subalignment_count = -1

    # Assuming FASTA format, collect all FASTA alignment files in the directory
    print('Collecting all pairwise alignments for concatenation!')
    for current_alignment in range(1, alignment_count):
        with open('{0}/seq_{1}.aln'.format(alignment_directory, str(current_alignment))) as fasta_input:
            # The subalignment_positions file is used to record the
            # start and stop positions of each aligned protein sequence in
            # the concatenation

            for line in fasta_input:
                stripped_line = line.strip()
                if '>' in stripped_line:
                    # Recognize the aligned species
                    if stripped_line[1:] == species1:
                        current_species = 0
                    elif stripped_line[1:] == species2:
                        current_species = 1
                else:
                    sequences[current_species] += stripped_line

            # Add the current sequence length to the subalignment positions
            # The protein ID cannot be added at this point, since protein IDs
            # are stripped to species IDs. You would need to record the
            # subalignment positions when reading in the sequences
            if subalignment_count > -1:
                subalignment_positions.append([subalignment_positions[subalignment_count][1] + 1,len(sequences[0]) - 1])
            else:
                subalignment_positions.append([0,len(sequences[0]) - 1])
            subalignment_count += 1

    # Perform a bootstrap analysis on the alignment and note the created variance
    def generate_bootstraps(alignment_length, generated_bootstrap_count):

        # Seed the singleton random module using the current system clock time
        # (by passing no parameter)
        seed = random.randrange(sys.maxsize)
        random.seed(seed)

        print('The seed for bootstrapping the alignment is: ', seed)
        with open('bootstrap_sample_rng_seed.R','w') as seed_file:
            seed_file.write('# The seed for sampling columns for pairwise distance\n# bootstrap analysis is a follows:\n# '+ str(seed) + '\n')

        # Generate a list of sequence indices to sample
        # The list is sampled from all alignment indices with equal weights and replacement
        for c in range(generated_bootstrap_count):
            yield random.choices(range(alignment_length), k=alignment_length)

#    if bootstrap_count > 0:
#        bootstrap_alignment_indices = create_bootstraps(len(sequences[0]),bootstrap_count)

    print('Concatenate and duplicate the pairwise aligned sequences!')
    # Generate the main concatenated and duplicated alignment for calculating
    # the pairwise species distance
    #sequence_pair_to_phylip(species1,species2,sequences,subalignment_positions)
    sequence_pair_to_phylip(species1,species2,sequences)

    def generate_multiprocessing_to_phylip_args_list(species1, species2, sequences, subalignment_positions, bootstrap_count):
        i = -1
        # Sample a set of position indices to draw aligned amino acids from both sequences
        for bootstrap_position_indices in generate_bootstraps(len(sequences[0]), bootstrap_count):
            i += 1
            yield [species1,species2,sequences,subalignment_positions,bootstrap_position_indices,i]

    if bootstrap_count > 0:
        print('Concatenate and duplicate the bootstrapped pairwise alignments!')
        # Generate the bootstrapped versions of the concatenated alignment for
        # estimating the distance variance
        # Each alignment bootstrap creates a separate process
        # The pool iterates through "bootstrap_alignment_indices"
        try:
            with Pool(nr_processors) as process_pool:
                # imap is a lazy loader of the bootstrap position indices
                # Sampling every position and putting their integers in a list takes
                # much more RAM than when the positions are written into two strings
                process_pool.imap(sequence_pair_to_phylip_multiprocessed, generate_multiprocessing_to_phylip_args_list(species1,species2,sequences,subalignment_positions,bootstrap_count),5)
        except KeyboardInterrupt:
            sys.exit('The user interrupted the generation of bootstrapped alignments!')

def search_protein_pairs(species1, species2, pair_table):
    """ Searches the given table for mutual presences of
    species 1 and species 2 within the first two columns
    and returns both columns verbatim in the order the
    species were passed to this function. """

    def generate_large_orthologous_pair_file_buffers(pair_table):
        # Read in oma_pair lines that contain proteins
        # of both species together
        with open(pair_table,'r') as pair_mapping:
            # The buffersize has been set that the OMA pairs file
            # approximately takes 1 GB per buffer on RAM
            tmp_lines = pair_mapping.readlines(134217728)
            while tmp_lines:
                # Identify the input species pair
                # among the orthologous protein pairs
                yield tmp_lines
                tmp_lines = pair_mapping.readlines(134217728)

    def generate_pairwise_orthologous_lines(species1, species2, pair_table):
        for buffered_lines in generate_large_orthologous_pair_file_buffers(pair_table):
            if buffered_lines:
                for line in buffered_lines:
                    if species1 in line and species2 in line:
                        yield line

    try:
        sequence_count = 1
        for line in generate_pairwise_orthologous_lines(species1, species2, pair_table):
            columns = line.rstrip().split('\t')
            # The previous generator already established that
            # both species are present in the columns. Here, we
            # look at the exact column position.
            if species1 in columns[0]:
                yield (columns[0], columns[1])
            else:
                yield (columns[1], columns[0])

            # Inform the user about the current progress
            if sequence_count % 1000 == 0:
                print('Sequence Nr.: {0}'.format(str(sequence_count)))
        print('Sequence Nr.: {0}'.format(str(sequence_count)))

    except KeyboardInterrupt:
        sys.exit('The user interrupted the compilation of orthologous pairs!')
    except FileNotFoundError:
        sys.exit('ERROR: The OMA pairs file is missing!')

def provide_protein_pair_sequences(species1, species2, sequence_source_file):
    """ Provide a dict of protein ids and sequences. """

    # All proteomes can be located in separate directories
    # that follow a systematic nomenclature
    # Otherwise, we assume that all proteins of the
    # species can be found within the oma_seqs.fa file
    #sequence_sources = set()
    #if os.path.exists(oma_proteomes_dir + '/proteome_' + species1):
    #    sequence_sources.add(oma_proteomes_dir + '/proteome_' + species1)
    #else:
    #    sequence_sources.add(oma_seqs)
    #if os.path.exists(oma_proteomes_dir + '/proteome_' + species2):
    #    sequence_sources.add(oma_proteomes_dir + '/proteome_' + species2)
    #else:
    #    sequence_sources.add(oma_seqs)

    try:
        sequence_index = SeqIO.index(sequence_source_file, 'fasta')
        if species1+"00001" not in sequence_index or species2+"00001" not in sequence_index:
            print('ERROR: Species are not represented in the sequence fasta file!')
            sys.exit()
        return sequence_index
    except KeyboardInterrupt:
        print('Gathering the sequences of pairwise orthologs was interrupted by the user.')
        sys.exit()
    except FileNotFoundError:
        print('ERROR: The FASTA file that contains the sequences is missing!')

def align_pairwise_proteins_pairwise2_parallelized(args):
    """ Provides a top-level function to parallelize the alignment of two sequences,
        args[0] and args[1]. Other parameters are static for every pair of sequences."""

    # one_alignment_only retrieves the first of the equally best scoring alignments.
    # d: match scores are derived from a list (the substitution matrix)
    # s: same gap penalties for both sequences (values for gapopen and gapextend were taken from MAFFT 6 L-INS-i
    # One alignment consists of a tuple of four elements:
    # sequence 1, sequence 2, score, beginning, 0-based character count

    return pairwise2.align.globalds(args[0], args[1], matlist.blosum62, -10.0, -1.0, penalize_end_gaps=True, one_alignment_only=True)[0]

def calculate_protein_distances_parallelized(args):
    """ A pickle-able dispatcher for calculate_protein_distances. """
    calculate_protein_distances(*args)

### Calculates the pairwise species maximum likelihood distance
### between species1 and species2. The distance is copied to the
### target_dir. The config is loaded from ProtTrace's configure.py
### if this script is executed directly.
def calculate_protein_distances(species1, species2, config, target_dir, preset_nr_processors = None, add_filename=1, bootstrap_count=0):

    print("Calculating the protein distance between {0} and {1}.".format(species1,species2))

    # Read the configuration of ProtTrace for paths
    oma_seqs = config.path_oma_seqs
    oma_pairs = config.path_oma_pairs
    concatAlignment = config.concat_alignments_script
    oma_proteomes_dir = config.path_distance_work_dir
    #linsi = config.msa
    treepuzzle = config.treepuzzle
    if treepuzzle is None:
        print('ERROR: No path to TreePUZZLE has been configured!')
        sys.exit()
    if preset_nr_processors is None:
        nr_processors = config.nr_processors
    else:
        nr_processors = preset_nr_processors
    delete_temp = config.delete_temp

    """ Directory management """

    root_dir = os.getcwd()
    # Create the working directory for the processed species pair
    work_dir = species1+'_'+species2
    if not os.path.exists(work_dir):
        os.mkdir(work_dir)
    # The directory which contains the fasta files
    # fa_dir = os.path.abspath(work_dir) + '/fa_dir'
    # The directory which contains the alignments
    # aln_dir = os.path.abspath(work_dir) + '/aln_dir'
    #if not os.path.exists(fa_dir):
    #    os.mkdir(fa_dir)
    #if not os.path.exists(aln_dir):
    #    os.mkdir(aln_dir)

    os.chdir(work_dir)

    # Align the set of pairwise oprthologs between the query
    # and the target species
    # print('Gather orthologous protein pairs.')
    #print('Preprocessing:\tParsing sequence pairs and aligning them...')


    # If the counter has never been incremented,
    # we can assume that the species pair is missing in the file
    #if sequence_count == 1:
    #    print("ERROR: No orthologous pairs found between {0} and {1}!".format(species1,species2))
    #    os.chdir(root_dir)
    #    return

    # print('Gather the sequences of all gathered orthologous protein pairs.')

    print('Aligning and concatenating pairwise orthologous protein sequences.')

    def align_pairwise_proteins_mafft(nr_processors, fa_dir, sequence_count):
        # We perform a global alignment of each orthologous pair
        # For this step, we only need to know the number of pairs
        # to access the fasta files in the fasta directory
        # ProtTrace rather uses MAFFT linsi than muscle

        # Encapsulate the call of MAFFT to enable parallelization
        def align_protein_pair(filename):
            # Execute the subprocess and wait for completion within the process
            # The wait function within a threaded function is necessary, because,
            # as far as I know, we can not manage a pool of threads with subprocess.Popen
            subprocess.Popen('{0} --amino --quiet --thread 1 {1} > {2}'.format(linsi, filename, filename.replace('.fa', '.aln').replace(fa_dir, aln_dir)),shell=True).wait()

        try:
            # Measure the time taken
            print('Aligning orthologous protein pairs')
            start = time.time()
            # To avoid adding another dependency, I uncommented the muscle command
            #os.system('muscle -quiet -in %s -out %s' %(filename, filename.replace('.fa', '.aln').replace(fa_dir, aln_dir)))
            # Run the alignment program in parallel threads
            tp = ThreadPool(nr_processors)
            tp.map(align_protein_pair,['{0}/seq_{1}.fa'.format(fa_dir,str(c)) for c in range(1,sequence_count)])
            tp.close()
            tp.join()
            # Print the time passed
            print('Total time passed for aligning orthologous proteins: ' + str(time.time() - start))

        except KeyboardInterrupt:
            print('The user has stopped the generation of alignments')
            tp.close()
            tp.join()
            sys.exit()

    def generate_pairwise2_alignment_arguments(prot_pairs_generator, sequence_records):
        """ Collects the sequences of pairwise orthologous proteins. """
        for prot_pair in prot_pairs_generator:
            yield (str(sequence_records[prot_pair[0]].seq), str(sequence_records[prot_pair[1]].seq))

    def align_pairwise_proteins_pairwise2(prot_pairs_generator, sequence_records, nr_processors):
        """ Aligns all pairwise orthologous protein sequences between species1 and species2 using
            the Biopython pairwise2 module. If the distance calculation between species is done
            sequentially, then this process can spawn child processes by providing a nr_processors
            count greater than 1. """

        if nr_processors > 1:
            try:
                with Pool(nr_processors) as alignment_process_pool:
                    # imap is a lazy loader of arguments. Very long sequences could take too much
                    # RAM. The unordered version is also a bit faster and does not need to preserve
                    # earlier, but later needed results in memory. With the last parameter, we allow
                    # imap to load in a couple of sequences.
                    for alignment in alignment_process_pool.imap_unordered(align_pairwise_proteins_pairwise2_parallelized, generate_pairwise2_alignment_arguments(prot_pairs_generator, sequence_records), 10):
                        yield alignment

            except KeyboardInterrupt:
                sys.exit('The user interrupted the alignment of pairwise orthologous proteins.')

        else:
            try:
                for arguments in generate_pairwise2_alignment_arguments(prot_pairs_generator, sequence_records):
                    # We use the parallelizable function here only to pass all arguments as a list.
                    yield align_pairwise_proteins_pairwise2_parallelized(arguments)

            except KeyboardInterrupt:
                sys.exit('The user interrupted the alignment of pairwise orthologous proteins.')
            # Since this code section will be executed if the species pairs are parallelized,
            # we need to explicity propagate any exception ot the parent process.
            except Exception as e:
                os.chdir(root_dir)
                raise e

    sequence_records = provide_protein_pair_sequences(species1, species2, oma_seqs)

    # The collection of orthologous pairs, their sequences, their alignments and their concatenation are all chained with generators here.
    concatenate_alignment(species1, species2, align_pairwise_proteins_pairwise2(search_protein_pairs(species1, species2, oma_pairs), sequence_records, nr_processors), bootstrap_count, nr_processors)

    #align_pairwise_proteins_mafft(nr_processors, fa_dir, sequence_count)

    # Collect all pairwise protein alignments, concatenate them
    # and calculate the summarized pairwise species distance
    #print('Preprocessing complete..\nConcatenating the alignments..')

    # Concatenate and bootstrap the pairwise protein alignments
    #concatenate_alignment_legacy(species1, species2, sequence_count, aln_dir, bootstrap_count, nr_processors)

    # The concatAlignment script requires perl
#    os.system('perl %s -in=%s -out=%s' %(concatAlignment, aln_dir, concat_file))
    #print('Postprocessing the concatenated alignment file..')

    # Executes TreePUZZLE to calculate the tree distance between the species pair
    def calculate_pairwise_distance(concat_file, treepuzzle):

        # Prepare the parameter file for TreePUZZLE
        with open('temp_puzzleParams.txt', 'w') as pp:
            pp.write(concat_file + '\ne\nm\nm\nm\nm\nm\nm\ny\n')

        # Execute TreePUZZLE
        os.system('{0} < temp_puzzleParams.txt >/dev/null'.format(treepuzzle))

    def read_calculated_distance(concatenated_protein_set_filename):
        try:
            with open(concatenated_protein_set_filename + '.dist','r') as concat:
                next(concat)
                return concat.readline().split()[3]
        except FileNotFoundError:
            sys.exit('Error: The pairwise distance calculation result file from TreePUZZLE cannot be found!')

    def write_main_result(output_filename, distance, target_dir, add_filename_to_target_dir):

        # Write the main computed pairwise species distance to a result file
        with open(output_filename,'w') as result:
            result.write(distance + '\n')

        # Copy the main output file to the given target directory, if given
        # The target directory is originally designed to be the ProtTrace cache directory
        if target_dir is not None:
            # If the target_dir is just a directory, append the file name
            if add_filename_to_target_dir == 1:
                target_dir_copy_path = os.path.join(os.path.abspath(target_dir), output_filename)
            # If the target_dir looks like a filename on its own, use it directly
            else:
                target_dir_copy_path = target_dir
            with open(target_dir_copy_path,'w') as result:
                result.write(distance + '\n')

    # Process the pairwise species maximum likelihood distance
    print('Calculating and writing the maximum likelihood distance!')
    concat_filename = species1 + '_' + species2 + '.phy'
    result_file = concat_filename.replace(".phy",".lik")
    calculate_pairwise_distance(concat_filename,treepuzzle)
    write_main_result(result_file,read_calculated_distance(concat_filename),target_dir,add_filename)

    # Process the distances calculated from bootstrapped alignments
    # First, we define a couple of generator functions

    # Anticipate the output file names of the
    # bootstrapped concatenated alignments
    def generate_bootstrap_concat_filenames(bootstrap_count):

        for i in range(0, bootstrap_count):
            yield species1 + '_' + species2 + '_bootstrap_' + str(i) + '.phy'

    def generate_bootstrap_distances(bootstrap_count, treepuzzle):

        for bootstrap_concat_filename in generate_bootstrap_concat_filenames(bootstrap_count):
            calculate_pairwise_distance(bootstrap_concat_filename,treepuzzle)
            yield read_calculated_distance(bootstrap_concat_filename)

    def calculate_and_write_bootstrap_distance_table_output(bootstrap_count, treepuzzle):

        # Write the main distance and its bootstrap distances into a table
        # for followup statistical analyses
        # The bootstrap column is a R friendly boolean to separate the main distance
        # from the bootstrap values
        separator = '\t'
        columns = [separator.join(['Species_1','Species_2','Distance'])]
        for distance in generate_bootstrap_distances(bootstrap_count,treepuzzle):
            columns.append(separator.join([species1,species2,distance]))
        with open(species1 + '_' + species2 + '_bootstrap.lik','w') as analysis_output:
            analysis_output.write('\n'.join(columns) + '\n')

    if bootstrap_count > 0:
        calculate_and_write_bootstrap_distance_table_output(bootstrap_count, treepuzzle)

    if delete_temp:
        delete_temporary_files(species1, species2, fa_dir, aln_dir, bootstrap_count, result_file)

    print('Finished species pair {0} - {1}'.format(species1,species2), flush = True)
    os.chdir(root_dir)

    # END OF POSTPROCESS FUNCTION DEFINITION

# This defines the start of this script if someone wants to calculate
# distances separately with this script
if __name__ == "__main__":

    # Define and check for arguments
    try:
        opts, args = getopt.getopt(argv, "q:s:o:r:h", ["query=", "species=", "output=", "help"])
    except getopt.GetoptError:
        print('Invalid arguments:\nUsage:\tdistanceCalculation.py -q <query_species_id> -s <paired_species_id> -o <output_dir> [-help]')
        sys.exit(2)

    for opt, arg in opts:
        if opt in ('-h','--help'):
            print("USAGE:\tdistanceCalculation.py -q <query_species_id> -s <paired_species_id> -o <output_file>\n\t-q\t\tThe OMA ID of the query species\n\t-s\t\tThe OMA ID of the target species\n\t-o\t\tThe full output file path")
            sys.exit(2)
        elif opt in ('-q', '--query'):
            query = arg
        elif opt in ('-s','--species'):
            species = arg
        elif opt in ('-o','--output'):
            output = arg
        else:
            print('Invalid arguments:\nUsage:\tprotTrace.py -i <omaIdsFile> | -f <fastaSeqsFile> -c <configFile> [-help]')
            sys.exit(2)

    if len(sys.argv[1:]) == 0:
        print('ERROR:\tNo arguments entered for calculating species distances:\nUSAGE:\tdistanceCalculation.py -s <query species ID> -c <configFile> [-help]')
        sys.exit(2)

    # If this script is started on its own calculate the pairwise distance between the specified species pair
    calculate_protein_distances(query,species,None,output)
