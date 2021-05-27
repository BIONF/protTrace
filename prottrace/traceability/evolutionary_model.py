#!/usr/bin/env python
# -*- coding: utf-8 -*-

#######################################################################
# Copyright (C) 2020 Arpit Jain,  Dominik Perisa,
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

# Preprocessing steps of ProtTrace include finding orthologs,  their sequences
# and computing their evolutionary model parameters.

import os
import sys

from subprocess import (
    run,
    CalledProcessError,
    PIPE
)
from pathlib import Path
from dendropy import TreeList

from traceability.tree_reconstruction import main as tree_reconst
from traceability.transform_alignment import main as trans_align
from utils.data_api import (
    orth_group
)
from utils.file import (
    dir_move
)
from utils.log import (
    print_progress,
    print_warning,
    print_error,
    time_report
)
from traceability.evol_params import evol_parameters


def main(query, prot_config):
    """ Reads in the proteome and calculates the evolutionary
    parameters of the given protein identifier. """

    # Print the current protein ID to the user screen
    print_progress(f'Currently processed protein ID: {query.id}')

    # The cache manager organizes calls to the ProtTrace cache dir.
    cache = cache_manager(prot_config)
    # Fetches the main species id from the configuration file
    # species_id = prot_config.species

    # Move to the current query's working directory
    work_dir = dir_move(prot_config.path_work_dir / query.id)

    # Instantialize the process timer and start the timer.
    timer = time_report()

    # Executing the steps of this module.
    exec_steps(query, prot_config, work_dir.current, cache)

    timer.print_time('Total preprocessing time.', scale='hrs')

    # Moving back to the previous working directory.
    work_dir.reverse()


def exec_steps(query, prot_config, work_dir_path, cache):
    """ Dispatches the steps of this preprocessing script. """

    nr_processors = prot_config.nr_processors

    """ Retrieve the orthologous sequences. """

    orth = orth_group(prot_config, query.id)

    # Retrieves orthology information from OMA database or infers
    # orthologs from hamstr / hamstr-OneSeq.
    exec_step(prot_config.orthologs_prediction,
              orth.filepath(query.id),
              "Orthologs file exist. Reusing it.",
              "Retrieving orthology information.",
              "Orthology information.",
              cache,
              orthology_prediction,
              query, prot_config)

    # Extra layer of security,  because the orthologs file must be
    # present for ProtTrace to run any further.
    if not orth.filepath(query.id).exists():
        print_error(('No orthlog file found in the working directory '
                     f'{str(work_dir_path)}. If the file is present in the '
                     'cache directory, please turn ON the cache flag in '
                     'program configuration file.'))
        sys.exit()

    """ Calculates evolutonary parameters. """
    """ The first step is to create a MSA of all orthologous sequences. """

    # Defines the filename of the MSA file in PHYLIP format.
    phy_file = Path(f'ogSeqs_{query.id}.phy').resolve()

    # Performs MSA on the orthologous sequences.
    exec_step(prot_config.perform_msa,
              phy_file,
              "Re-using previously compiled orthologs alignment.",
              "Performing MSA of the orthologous sequences",
              "MAFFT",
              None,
              multiple_sequence_alignment,
              prot_config.msa, nr_processors, orth.path, phy_file, 'phylip')

    """ Initialize the calculation of evolutionary parameters. """

    # The constructor inserts default values from the configuration.
    evol_params = evol_parameters(query, prot_config)

    """ Calculation of the substitution rate scaling value. """

    # Defines the file name of the reconstructed species tree. This tree
    # is used to extract the mean phylogenetic speed between species.
    # But also later when normalizing indel rates.
    tree_file = Path(f'ogSeqs_{query.id}.phy.treefile')

    # Reconstructs the phylogeny of the orthologous sequences and
    # calculates the substitution rate scaling from this tree.
    exec_step(prot_config.calculate_scaling_factor,
              None,
              "Pre-computed scaling factor found for re-use!",
              "Tree reconstruction and scaling factor calculation.",
              "RAxML",
              cache,
              tree_reconst,
              query, orth, evol_params, tree_file, cache, prot_config)

    """ Calculation of the indel rate and indel length distribution value. """

    # Calculates indel evolutionary values. Depending on the ProtTrace
    # configuration,  values are computed assuming parsimony or using
    # SpartaABC.
    exec_step(prot_config.calculate_indel,
              None,
              'Pre-computed indel values found for re-use!',
              'Indel rate and length distribution calculation.',
              'Indels',
              cache,
              calculate_indels,
              query, evol_params, tree_file, phy_file, orth, prot_config)

    """ Compiles preprocessing steps into the REvolver configuration file. """

    # Checks necessary file paths before expensive functions.
    if not prot_config.pfam_database.exists():
        print_error(f'The file {str(prot_config.pfam_database)} does not '
                    'exist. Check the path to the Pfam database file.')
        sys.exit()
    if not prot_config.hmmfetch.exists():
        print_error(f'The file {str(prot_config.hmmfetch)} does not exist. '
                    'Check the path to the hmmfetch program.')
        sys.exit()

    # Here we create the domain constraint file for REvolver.

    # The HMM file of the orthologous proteins.
    hmm_file = query.id + '.hmm'

    # Execution of hmmscan.
    print_progress("Generating domain constraints for REvolver!")
    run_hmmscan(query, hmm_file, prot_config)

    # Prepare the XML configuration file for REvolver.
    prepare_REvolver_XML(prot_config, query, evol_params, work_dir_path)


def prepare_REvolver_XML(config, query, evol_params, work_dir_path):
    """ Writes the configuration file for REvolver. """

    xml_file = Path(f'revolver_config_{query.id}.xml')

    pfam_db = config.pfam_database
    hmmfetch = config.hmmfetch
    aa_matrix = config.aa_substitution_matrix
    indel_rate = evol_params.indel_rate
    indel_dist_p = evol_params.indel_length_distribution
    scaling_factor = evol_params.scaling_factor
    evol_tree = config.simulation_tree
    hmm_file = f'{query.id}.hmm'
    output_dir = work_dir_path / 'REvolver_output'
    if not output_dir.exists():
        output_dir.mkdir()
    run_spartaABC = config.run_spartaABC

    if float(indel_dist_p) > 1.0:
        if not run_spartaABC:
            print_warning('The indel length distribution parameter is greater '
                          'than 1, indicating a zipfian slope used by '
                          'SpartaABC. This contradicts the run_spartaABC '
                          'option being turned off')
        length_dist_section = ('\t\t\t\t<length distribution="zipfian_M" '
                               f'z="{indel_dist_p}" max="50"/>\n')
    else:
        length_dist_section = ('\t\t\t\t<length distribution="geometric" '
                               f'p="{indel_dist_p}"/>\n')

    # We produce the XML file in a litteral string. This way, we do not need
    # to utilize a libary together with the risk of version dependent changes.
    xml_c = (''
             '<configdata  xsi:schemaLocation="http://www.cibiv.at/Revolver '
             './input_schema.xsd" xmlns="http://www.cibiv.at/Revolver" '
             'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" >\n'
             '\t<config>\n'
             f'\t\t<hmmdb path="{pfam_db}"/>\n'
             f'\t\t<hmmfetch location="{hmmfetch}"/>\n'
             '\t</config>\n'
             '\t<model>\n'
             f'\t\t<substitution name="{aa_matrix}"/>\n'
             '\t\t<indel>\n'
             f'\t\t\t<insertion rate="{indel_rate}">\n'
             f'{length_dist_section}'
             '\t\t\t</insertion>\n'
             f'\t\t\t<deletion rate="{indel_rate}">\n'
             f'{length_dist_section}'
             '\t\t\t</deletion>\n'
             '\t\t</indel>\n'
             '\t</model>\n'
             f'\t<tree scalingFactor="{scaling_factor}" '
             f'path="{str(evol_tree)}"/>\n'
             '\t<root>\n'
             '\t\t<inputSequence>\n'
             f'\t\t\t<fasta file="seq_{query.id}.fa"/>\n'
             f'\t\t\t<hmmer file="{hmm_file}"/>\n'
             '\t\t</inputSequence>\n'
             '\t</root>\n'
             '\t<output>\n'
             f'\t\t<dir path="{output_dir}" separateFastaFiles="false" '
             'trueAlignment="false" include="leaf"/>\n'
             '\t</output>\n'
             '</configdata>'
             )

    with xml_file.open('w') as xf:
        xf.write(xml_c)


def run_hmmscan(query, hmm_file, config):
    """ Writes the query sequence to a fasta file. Hmmscan will then create
    a pHMM domain constraint file out of it for REvolver by detecting Pfam
    domains. """

    # This is also important for REvolver later!
    seq_file = query.write_to_fasta()

    hmmscan = config.hmmscan
    pfam_database = config.pfam_database
    hmm_file = f'{query.id}.hmm'

    # Execute hmmscan.
    try:
        run([hmmscan, '--notextw', '-o', hmm_file, '-E', '0.01',
             pfam_database, seq_file], check=True)
    except CalledProcessError:
        print_error('Hmmscan did not run correctly. Check the given pfam '
                    'files and the query sequence')
        sys.exit()


def orthology_prediction(query, prot_config):
    """ Predicts orthologs to the query sequence among the species
    in the species mapper. """

    try:
        return orth_group(prot_config, query.id)
    except KeyboardInterrupt:
        sys.exit(print_error("Keyboard interruption by user!"))

# Calculate evolutionary parameters
# Substitution rate scaling is outsourced to treeReconstruction.py
# Therefore,  most stuff here is for calculating the indel rates and
# indel length distribution.


def multiple_sequence_alignment(program, nr_processors, orth_file, output_file,
                                output_format):
    """ Performs a multiple sequence alignment of the protein sequences
    within the orthologs file. This function accepts fasta and phylip
    as output formats. """

    try:
        # The command arguments for MAFFT L-INS-I in subprocess format.
        # Element 2 is reserved for specifying the output format (PHYLIP/FASTA)
        cmd = [str(program), '--quiet', '', '--thread', str(nr_processors),
               str(orth_file), '>', str(output_file)]

        # If the output format parameter explicitly demands the MSA
        # to be outputted into PHYLIP format, other than FASTA,
        # the necessary format argument is added to the subprocess call.
        if output_format == 'phylip':
            cmd[2] = '--phylipout'
        else:
            del cmd[2]

        # Runs the MSA performing program. The pipe for MAFFT requires
        # shell=TRUE, I would like to be proven wrong. Shell=TRUE also requires
        # the command to be passed as a single string.
        run(' '.join(cmd), check=True, shell=True)

    except CalledProcessError:
        # Failed MSAs result in a warning.
        print_warning('MSA did not work. Probably less than 2 sequences were '
                      'found for this alignment')


def calculate_indels(query, evol_params, tree_file, phy_file, orth,
                     prot_config):
    """ Calculate the indel rate and indel length distribution.
        This master function processes the configuration to compute
        them either with SpartaABC or assuming parsimony. """

    # Calculating indel values require a precomputed protein tree in
    # any case. If no tree is found, we still use default indel values.
    if tree_file.exists():
        # Decide whether SpartaABC is used for calculating indel values.
        if prot_config.run_spartaABC:
            path_posterior_params = Path(
                'spartaABC_raw_output.posterior_params')

            # Calculate indel values using SpartaABC.
            exec_step(prot_config.run_sparta_abc,
                      path_posterior_params,
                      'Re-using already computed posterior estimates.',
                      'Calculating indel parameters with SpartaABC',
                      'SpartaABC',
                      None,
                      calculate_indels_sparta_abc,
                      query, evol_params, orth, tree_file, prot_config)
        # The alternative to SpartaABC is to calculate indel values
        # assuming parsimony. A parsimonious gap model is parsed to
        # IQTree --tina.
        else:
            calculate_indels_parsimony(query, evol_params, phy_file, tree_file,
                                       prot_config)


def calculate_indels_sparta_abc(query, evol_params, orth, tree_file,
                                prot_config):
    """ Calculates indel rate and indel length distribution using SpartaABC.
    """

    # Create a new MSA of the orthologous sequences,  but in
    # FASTA format.

    # Define the filename for the MSA file.
    aln_file = Path(f'ogSeqs_{query.id}.fa')

    # Perform the MSA using MAFFT.
    exec_step(prot_config.perform_msa,
              aln_file,
              "Re-using previously compiled orthologs alignment.",
              "Performing MSA of the orthologous sequences",
              "MAFFT",
              None,
              multiple_sequence_alignment,
              prot_config.msa, prot_config.nr_processors, orth.path,
              aln_file, 'fasta')

    # Parse important parts of the ProtTrace configuration.
    path_spartaabc = prot_config.sparta
    # DAWG performs on nucleotides, which is unsupported by ProtTrace.
    # DAWG support is considered deprecated.
    # seq_evolve_with_dawg = prot_config.evolve_dawg
    path_posterior_params = Path('spartaABC_raw_output.posterior_params')

    # The delete temporary files option needs both control files
    # to be named to check them.
    indelible_control_file = "indelible_control.txt"
    dawg_control_file = "dawg_control.txt"

    # Implementation of SpartaABC's helper script to determine
    # little less than min sequence length and little more than
    # the max sequence length
    RL_values = calculate_sparta_rl_values(orth)

    # Both sequence evolution programs are to be specified inside
    # the config file. The left out one is specified with blanks.
    indelible_control_file = ""
    dawg_control_file = ""

    # The entire newick tree of the orthologous group must be
    # written inside the respective control file.
    written_out_tree = ""
    # Try to read in the precomputed orthologous sequence tree.
    try:
        with open(tree_file, 'r') as tree:
            written_out_tree = tree.read().rstrip('\n')
    except IOError:
        sys.exit(
            print_error('There is no orthology tree file available! Make sure '
                        'the option orthologs_tree_reconstruction has been '
                        'executed before.'))

    # Write the necessary control file.
    # if seq_evolve_with_dawg:
    #     # Uses Dawg for sequence evolution.
    #     dawg_simulator = '1'
    #     # The INDELible control file specification is left empty.
    #     dawg_control_file = 'dawg_control.txt'
    #     # A newline character is added manually at the end.
    #     control_file_contents = '\n'.join((
    #         '[Tree]',
    #         'Tree = \'{0}\''.format(written_out_tree),
    #         '',
    #         '[Indel]',
    #         'Subst.Model = \'WAG\'',
    #         'Model.Ins = \'POWER-LAW\'',
    #         'Model.Del = \'POWER-LAW\'',
    #         'Rate.Ins = ?',
    #         'Rate.Del = ?',
    #         'Params.Ins = ?,  50.0',
    #         'Params.Del = ?,  50.0',
    #         'Max.Ins = 50.0 # necessary to avoid bad_alloc error for now',
    #         'Max.Del = 50.0',
    #         '',
    #         '[Sim]',
    #         'Reps = 1',
    #         '[Output]',
    #         'File = \'fasta:-\'',
    #         '[Root]',
    #         'Length = ?',
    #         ''))

    #     fill_file(dawg_control_file,  control_file_contents)

    else:
        # Uses INDELible for sequence evolution.
        dawg_simulator = '0'
        # The Dawg control file specification is left empty
        indelible_control_file = 'indelible_control.txt'
        # A newline character is added manually at the end.
        control_file_contents = '\n'.join((
            '//  INDELible control file',
            '',
            '[TYPE] AMINOACID 2',
            '',
            '[SETTINGS]',
            '    [output] FASTA',
            '    [fileperrep] TRUE',
            '',
            '[MODEL]    modelname',
            '  [submodel]    WAG',
            '  [indelmodel]  POW  ? 50',
            '  [indelrate]   ?',
            '',
            '[TREE] treename {0}'.format(written_out_tree),
            '',
            '[PARTITIONS]   partitionname',
            '    [treename modelname  1]',
            '',
            '[EVOLVE]     partitionname  1  ?',
            ''))

        fill_file(indelible_control_file,  control_file_contents)

    # Writes the SpartaABC config file, using SpartaABC's amino acid
    # specific values.
    # A newline character is added manually at the end.
    config_file_contents = (
        f'_indelibleTemplateControlFile {indelible_control_file}\n'
        f'_dawgTemplateControlFile {dawg_control_file}\n'
        f'_dawgSimulator {dawg_simulator}\n'
        f'_inputRealMSAFile {aln_file}\n'
        f'_outputGoodParamsFile {str(path_posterior_params)}\n'
        '_numberOfSamplesToKeep 100000\n'
        '_alignmentMode 0\n'
        '_similarity_mode 0\n'
        f'_minRLVal {RL_values[0]:.1f}\n'
        f'_maxRLVal {RL_values[1]:.1f}\n'
        '_minIRVal 0.0\n'
        '_maxIRVal 0.15\n'
        '\n'
        '_wAvgUniqueGapSize 0.0608829560660657\n'
        '_wMSAMin 0.00191213860761342\n'
        '_wNumGapsLenTwo 0.000864182910616497\n'
        '_wAvgGapSize 0.161663134960216\n'
        '_wTotNumGaps 0.000139867523540041\n'
        '_wNumGapsLenAtLeastFour 0.000273522492154681\n'
        '_wNumGapsLenOne 0.000439656531427027\n'
        '_wMSAMax 0.00160268570388424\n'
        '_wMSALen 0.000241093525341335\n'
        '_wTotNumUniqueGaps 0.00115596033397199\n'
        '_wNumGapsLenThree 0.00142589360646065\n')

    fill_file('Sparta.conf',  config_file_contents)

    # Runs SpartaABC
    try:
        run([str(path_spartaabc), 'Sparta.conf'], check=True)
    except CalledProcessError:
        path_posterior_params.unlink()
        sys.exit(print_error('SpartaABC threw an error! Removing the '
                             'incomplete posterior parameters file.'))

    # Whether a new instance of SpartaABC was run or not,  extract the
    # posterior parameters. The function is a copy of the script that
    # is provided with the SpartaABC installation.
    posterior_estimates = extract50closest_distance_posterior_mean(
        path_posterior_params)
    indel = posterior_estimates[0]
    indel_distribution = posterior_estimates[1]

    # The indel length distribution is intended to describe a Zipfian
    # distribution. This distribution
    if indel_distribution < 1:
        print_progress('Indel distribution estimated to be lower than 1. '
                       'This is not compatible to REvolver Zipfian '
                       'distribution. Default of 1.1 is used.')
        indel_distribution = 1.1

    evol_params.set_indel_rate(indel)
    evol_params.set_indel_distribution(indel_distribution)

    # Cleanup temporary files. Errors in this section should not stop
    # the entire program.
    if prot_config.delete_temp:
        try:
            if path_posterior_params.exist():
                path_posterior_params.unlink()
            os.system('rm Sparta.conf')
            if (indelible_control_file != ""
                    and os.path.exists(indelible_control_file)):
                os.system('rm indelible_control.txt')
            if dawg_control_file != "" and os.path.exists(dawg_control_file):
                os.system('rm dawg_control.txt')
        except Exception:
            print_warning('Deletion of temporary files after SpartaABC '
                          'failed. ProtTrace continues.')


def calculate_sparta_rl_values(orth):
    """ Implementation of the script file distributed with SpartaABC
    to calculate the less than minimum and maximum sequence lengths. """

    min_length = 0
    max_length = 0
    with orth.path.open('r') as orthologs_list:
        for count, line in enumerate(orthologs_list, start=1):
            if count == 1 and line[0] != ">":
                raise IOError('Orthologs Fasta file does not start with a '
                              f'Fasta header! Affected file: {str(orth.path)}')
            if count % 2 == 0:
                if count == 2:
                    min_length = len(line)
                    max_length = len(line)
                else:
                    if len(line) > max_length:
                        max_length = len(line)
                    if len(line) < min_length:
                        min_length = len(line)

    less_than_min = int(min_length - 0.1*min_length)
    if less_than_min < 10:
        less_than_min = 10
    more_than_max = int(max_length + 0.1*max_length)

    return (less_than_min, more_than_max)


def extract50closest_distance_posterior_mean(sparta_output_file):
    """ After simulating the evolution from the MSA and the phylogenetic
    tree,  the 50 sequences with the closest distance to the input are
    considered for extracting the indel rates and indel length
    distributions. From these 50 sequences the mean values are drawn. """

    simulations = []
    with sparta_output_file.open('r') as raw:
        # First line are headers,  second line represents initial MSA
        simulations = raw.readlines()[2:]
    # Extract the first (distance) and fourth (indel rate) element of each
    # tab-delimited line. Then sort by distance (element 0 of resulting pair
    # list).
    distance_indel_zipfian_slope = [[float(line.split('\t', 1)[0]),
                                     float(line.split('\t', 4)[3]),
                                     float(line.split('\t', 4)[2])]
                                    for line in simulations]
    distance_indel_zipfian_slope.sort(key=lambda x: x[0])

    # Sum up the 50 lowest distance's indel rates
    fifty_best_posteriors = [i[1:] for i in distance_indel_zipfian_slope[:50]]
    fifty_best_indel_sum = 0
    fifty_best_zipfian_slope_sum = 0
    for distance_sorted_indel_zipfian_slope in fifty_best_posteriors:
        fifty_best_indel_sum += distance_sorted_indel_zipfian_slope[0]
        fifty_best_zipfian_slope_sum += distance_sorted_indel_zipfian_slope[1]

    # Return the mean of the 50 closest distance's indel rates
    return (fifty_best_indel_sum / 50, fifty_best_zipfian_slope_sum / 50)


def calculate_indels_parsimony(query, evol_params, phy_file, tree_file,
                               prot_config):
    """ Calculates indel rates under the assumption of parsimony
    with IQTree -tina. """

    """ Prepare the alignment to be read by IQTree. """

    print_progress('Transforming MSA based on indel blocks')

    # Try to transform the multiple sequence alignment. If this
    # throws an error,  continue with the default value.
    trans_file = Path(f'ogSeqs_{query.id}.trans')
    alignment_length = 0
    try:
        alignment_length = trans_align(phy_file, trans_file)
    except Exception:
        print_warning('An error has occurred when transforming the '
                      'alignment for parsimony based evolutionary '
                      'parameter computation')

    """ Read the tree file for calculating indel rates. """

    print_progress('Calculating indels based on parsimony.')

    # Reads the contents of the protein tree that was reconstructed for
    # calculating substitution scaling rates.
    try:
        trees = TreeList.get_from_path(tree_file, 'newick')
        tree_lengths = [tree.length() for tree in trees]
    except Exception as e:
        raise e
        print_error('Could not read the reconstructed species tree')
        sys.exit()

    """ Calculate the indel rate and indel length distribution of the query
    protein. """

    result = ''

    # Estimates indel rate and indel length distribution.
    try:
        command = [str(prot_config.iqtree), "-s", str(trans_file),
                   str(tree_file), "-tina", "-st", "MULTI"]
        command_str = ' '.join(command)
        print_progress(f'IQ-Tree command: {command_str}')
        result = run(command, stdout=PIPE, check=True).stdout.decode('utf-8')
    except CalledProcessError:
        print_warning('IQ-Tree could not estimate indel rates!')

    # Parses the result from IQ-Tree.
    for line in result.split('\n'):
        if line.split(':')[0] == 'mean length':
            # If the file lists a new mean indel length,  the listed
            # value replaces the default indel length distribution.
            if float(line.split(':')[1].replace(' ',  '')) > 0:
                indel_distribution = 1 / float(line.split(':')[1]
                                               .replace(' ',  ''))
                if indel_distribution >= 1:
                    indel_distribution = 0.99
                elif indel_distribution < 0.02:
                    indel_distribution = 0.02
        # The listed indel rate replaces the default indel rate.
        elif line.split(':')[0] == 'Parsimony score is':
            indel = (float(line.split(':')[1].replace(' ',  '')) /
                     (alignment_length * tree_lengths[0])) / 2

    # Informs the user about the estimated indel parameters.
    print_progress(f'Indel rate: {indel}; Indel length distribution: '
                   f'{indel_distribution}')

    evol_params.set_indel_rate(indel)
    evol_params.set_indel_length_distribution(indel_distribution)


""" Parse the query proteome, orthologous proteins and their sequences from
    OMA. """

# REWRITTEN INTO THE PROTEOME OBJECT IN DATA_API.py
# def parse_proteome(query, species_id,  proteome_file,  cache,  prot_config):
#     """ Parses the proteome of the specified species identifier either
#     from the OMA seqs file or the hamstr genome directory,  if it does
#     not exist inside the cache directory yet. Then,  create a BLAST
#     database from it for future BLAST searches. """
#
#     # Initializing config parameters
#     path_oma_seqs = prot_config.path_oma_seqs
#     search_oma_database = prot_config.search_oma_database
#     path_makeblastdb = prot_config.makeblastdb
#     hamstr = hamstr_API(prot_config)
#     path_hamstr_mapping = prot_config.species_map
#
#     # Checks if the proteome has already been parsed and saved inside
#     # the cache directory,  if the cache directory shall be reused.
#     if cache.in_cache(proteome_file):
#         print_progress("Proteome exist. Reusing it.")
#     else:
#         if search_oma_database:
#             protome_datasource = prot_config.path_oma_seqs
#         else:
#
#             proteome_datasource = hamstr.GetGenomePath(hamstr_id)
#         if search_oma_database:
#             # Searches the oma seqs file for the proteome sequences.
#
#             print_progress(
#                 "Parsing gene set for species {0} from OMA database".format(species_id))
#             species_found_flag = False
#
#             # The future contents of the proteome file are filled into
#             # this variable before writing the file.
#             prot_file_content = []
#             with open(path_oma_seqs, 'r') as seqs_file:
#                 for line in seqs_file:
#                     # If the species identifier is found again inside the OMA-Seqs file,
#                     # set the sequence found flag to true and write the header
#                     # and the sequence in the next line into the proteome file.
#                     if line[:6] == ">" + species_id:
#                         species_found_flag = True
#                         prot_file_content.append(line + next(seqs_file))
#
#             # Fills the proteome file with the sequences of the query
#             # protein. We did not trim the contents of any newline
#             # characters. Therefore,  the contents are joined together
#             # without any delimiter. The prot_file_content variable
#             # becomes deleted.
#             fill_file(proteome_file,  "".join(prot_file_content))
#
#             # Throw an error,  if the species was not found within the OMA database.
#             if not species_found_flag:
#                 sys.exit(print_error("Species {0} was not found inside the oma sequence file. Please make sure that the species exists within the OMA database. If not,  turn off the search_oma_database option in the prog.config file.".format(species_id)))
#         else:
#             # Searches the hamstr genome directory for the proteome sequences.
#
#             print_progress("Parsing gene set for species {0} from local (hamstr) BLAST directory".format(species_id))
#             # If the species OMA id is found in the last column of any line inside the species_tree_mapping file,
#             # it's corresponding first column must contain the hamstr ID
#             for splitted_line in generate_splitted_lines(path_hamstr_mapping):
#                 if species_id == splitted_line[-1]:
#                     hamstr_id = splitted_line[0]
#                     break
#
#             # Specify the genome directory in which the proteome can be found in.
#             genome_dir = hamstr.ResolvePath("genome_dir")
#             # Build the actual proteome file directory
#             species_genome_file = "{0}/{1}/{1}.fa".format(genome_dir, hamstr_id)
#             # If the proteome file exists inside the genome directory of hamstr,  copy the contents to the current directory.
#             if os.path.exists(species_genome_file):
#                 os.system('cp {0} {1}'.format(species_genome_file, proteome_file))
#             else:
#                 sys.exit(print_error("Reference species {0} was not found in the local (hamstr) BLAST directory!".format(species_id)))
#
#         # Copy the written proteome file into the cache directory.
#         cache.copy_to_cache(proteome_file)
#         prepare_blastdb(config,  proteome_file)
#
#
# def prepare_blastdb(config,  proteome_file):
#     """ Creating a BLAST database for future searches. """
#     print_progress('Making a BLAST database of the gene set.')
#     os.system('{0} -in {1} -input_type fasta -dbtype prot'
#               .format(config.path_makeblastdb,  proteome_file))


# REWRITTEN INTO THE ORTH_GROUP OBJECT IN DATA_API.PY
# def find_oMAGroup(query, query_seq, work_dir_path, proteome_file, prot_config):
#     """ Reads the OMA orthologous groups file and parses the ortholog
#     list for input OMA id """
#
#     # Definition of some variables that would bloat the parameter list.
#     id_file = "ogIds_{0}.txt".format(query.id)
#     species_id = prot_config.species
#     path_oma_group = prot_config.path_oma_group
#     path_oma_seqs = prot_config.path_oma_seqs
#     path_makeblastdb = prot_config.makeblastdb
#
#     # This clause determines,  whether only an OMA ID without a sequence
#     # is given,  or a fasta file with a sequence,  but not with a
#     # guaranteed OMA ID in it's header.
#     # In this case,  a query protein ID is given with no sequence.
#     if query_seq == 'None':
#         print_progress("OMA IDs given")
#         # OMA ID present == 2
#         id_group_annotation_status = 2
#         print_progress("Searching OMA ortholog group for given OMA id {0}".format(query.id))
#         # Write the OMA ID of the query protein into the omaID file.
#         with open (work_dir_path + "/omaID.txt", 'w') as oma_ID_File:
#             oma_ID_File.write(query.id)
#
#         id_file_content = ""
#         # Open the OMA ID file where all orthologous OMA IDs are stored in.
#         written = False
#         try:
#             # Iterate the OMA-groups.txt file,  where all OMA
#             # orthologous groups are stored.
#             for line in generate_lines(path_oma_group):
#                 # Hashtags mark the header of an OMA file.
#                 if not "#" in line and query.id in line:
#                     # Orthologous OMA IDs are listed in wide tabular
#                     # format. The orthologous sequence IDs are located
#                     # behind the first two columns. Join them with
#                     # newlines.
#                     id_file_content = '\n'.join(split_line(line)[2:])
#                     # We assume only one orthologous group to be
#                     # assigned to this protein. Therefore,  we stop after
#                     # the first hit.
#                     written = True
#                     break
#             # If no orthologous groups has been found. The query
#             # protein ID is considered the only ortholog.
#             if not written:
#                 id_file_content = query.id
#         except IOError:
#             sys.exit(print_error("Cannot find the OMA orthologous",
#                 "groups file. Check OMA files given as input!"))
#
#         # Fills the ortholog ID file with orthologous sequences' IDs.
#         fill_file(id_file,  id_file_content)
#     # In this case,  a query protein is given in FASTA format with their
#     # sequence.
#     else:
#         # OMA ID must be estimated == 1
#         id_group_annotation_status = 1
#
#         # We search for the OMA ID by sequence similarity within the
#         # OMA-seqs.fa file.
#         sequence_found_flag = False
#         try:
#             # Open the OMA-seqs.fa file.
#             with open(path_oma_seqs) as f:
#                 for line in f:
#                     # Testing if the species ID is found in the
#                     # sequence header.
#                     if ">" in line and species_id in line:
#                         # Assuming the OMA-seqs.fa file has been gotten
#                         # rid of in-sequence newlines.
#                         if TrimLine(f.next()).replace("*",  "") == TrimLine(query_seq).replace("*",  ""):
#                             # Found the identical sequence
#                             query.id_temp = trim_line(line)[1:]
#                             sequence_found_flag = True
#                             break
#         except IOError:
#             sys.exit(print_error('Cannot find the OMA sequence file. '
#                 'Check OMA files given as input!'))
#         # If the exact same sequence could not be found in the OMA
#         # sequence file,  a BLAST search is attempted.
#         if not sequence_found_flag:
#             print_warning('The input sequence is not present in '
#                 'the OMA sequence file. Performing a BLAST search to get '
#                 'the closest best hit.')
#             # Create a temporary directory for the BLAST search.
#             temp_dir = "temp_blast_{0}".format(query.id)
#             if not os.path.exists(temp_dir):
#                 os.mkdir(tmp_dir)
#             # Changes into the temporary directory.
#             os.chdir(tmp_dir)
#             # Copies the proteome into the temporary directory.
#             os.system('ln -s {0} proteome.fa'.format(proteome_file))
#             # Casts the proteome into a BLAST database to search against.
#             os.system('{0} -in proteome.fa -dbtype prot'.format(path_makeblastdb))
#             # Prepare an input file for BLAST containing the query sequence.
#             with open("temp_input_seq.fa", 'w') as temp_input:
#                 temp_input.write(">{0}\n{1}".format(query.id, query_seq))
#             # Performs the BLASTP search
#             subprocess.call([path_blastp,    '-query',  'temp_input_seq.fa',
#                 '-db',   'proteome.fa',  '-evalue',  '0.00001',  '-outfmt',  '6',
#                 '-out',   'temp.txt'])
#             #LEGACY-code
#             #os.system('{0} -query temp_input_seq.fa -db proteome.fa -evalue 0.00001 -outfmt 6 -out temp.txt'.format(path_blastp))
#             if os.path.exists("temp.txt"):
#                 with open("temp.txt", 'r') as blastp_output:
#                     if len(blastp_output.read().split("\n")) > 1:
#                         # The identifier in the first column of the first line is taken as the protein identifier.
#                         query.id_temp = split_line(blastp_output.read())[1]
#                         sequence_found_flag = True
#                     else:
#                         sys.exit(print_error("The BLASTP output file does not have more than 1 line! Check the output file at {0}/{1}/temp.txt!".format(query.id, temp_dir)))
#             else:
#                 sys.exit(print_error("The BLASTP search did not yield an output file!"))
#             # Change back to the general work directory of this protein.
#             os.chdir(work_dir_path)
#             # Delete the temporary directory,  if the option has been set.
#             if prot_config.delete_temp:
#                 os.system('rm -rf {0}'.format(temp_dir))
#
#         # If the sequence has been found either by sequence identity or
#         # a BLASTP search.
#         if sequence_found_flag:
#             id_group_annotation_status = 2
#             # Write the acquired OMA ID into the OMA ID file.
#             with open (work_dir_path + "/omaID.txt", 'w') as oma_ID_File:
#                 oma_ID_File.write(query.id_temp)
#
#             print_progress("Searching OMA ortholog group for {0}".format(query.id))
#             id_file_content = ""
#             # Open the OMA ID file where all orthologous OMA IDs are
#             # stored in.
#             with open(id_file, 'w') as id_file_file:
#                 written = False
#                 try:
#                     # Open the OMA-groups.txt file,  where all OMA
#                     # orthologous groups are stored.
#                     for line in generate_lines(path_oma_group):
#                         # Hashtags mark the header of an OMA file.
#                         if not "#" in line and query.id_temp in line:
#                             # Orthologous OMA IDs are listed in wide tabular
#                             # format. The orthologous sequence IDs are located
#                             # behind the first two columns. Join them with
#                             # newlines.
#                             id_file_content = '\n'.join(split_line(line)[2:])
#                             # We assume only one orthologous group to be
#                             # assigned to this protein. Therefore,  we stop
#                             # after the first hit.
#                             written = True
#                             break
#                     # If no orthologous groups has been found.
#                     if not written:
#                         id_group_annotation_status = 1
#                         id_file_content = query.id_temp
#                 except IOError:
#                     sys.exit(print_error("Cannot find the OMA orthologous",
#                     " groups file. Check OMA files given as input!"))
#
#                 # Fills the ortholog ID file with orthologous sequences' IDs.
#                 fill_file(id_file,  id_file_content)
#
#     # Return the ending status of this function.
#     # 1: OMA ID has been found,  but no sequence or OMA group.
#     # 2: OMA ID has been found together with a sequence. Indicates a
#     # successful run.
#     return id_group_annotation_status


# REWRITTEN INTO THE PIPELINE PROTEOME -> PROT_ID FOR THE QUERY AND EACH
# ORTHOLOGOUS GROUP MEMBER
# def find_oma_sequences(query,  prot_config):
#     """ Searches sequences for all OMA IDs in the orthologous group. """
#
#     # Initializing config parameters
#     path_oma_seqs = prot_config.path_oma_seqs
#     path_hamstr_mapping = prot_config.species_map
#
#     try:
#         # Create a pass filter for all species in the mapping file. We
#         # expect the OMA ID of each species to be the last column of every
#         # line.
#         species_id_filter = [splitted_line[-1]
#             for splitted_line in generate_splitted_lines(path_hamstr_mapping)]
#     except IOError:
#         sys.exit(print_error('Cannot find the species mapping file. '
#         'The default name is speciesTreeMapping.txt.'))
#
#     print_progress(f'Searching OMA ortholog sequences for {query.id}')
#     # Predefine the content list of the ortholog sequence file.
#     orth_file_content = []
#     # Read in all ortholog IDs
#     with open(id_file,  'r') as id_input:
#         orthologs_ids = id_input.read().split("\n")
#     # This counter keeps track of how many orthologous sequences have
#     # been found yet. If all IDs are found. The search ends.
#     ids_count = len(ids)
#
#     try:
#         # Iterate through the oma_seqs.fa file.
#         with open(path_oma_seqs,  'r') as sequence_file:
#             for line in sequence_file:
#                 # Check whether the line contains a sequence ID and
#                 # whether the id is in the list of orthologous IDs and
#                 # the mapping file.
#                 if (line[0] == '>'
#                         and line.split('\n')[0][1:] in orthologs_ids
#                         and line[1:6] in species_id_filter):
#                     # Write the ID and sequence down. Remove undefined
#                     # amino acids. These are coded with either '*' or 'X'.
#                     orth_file_content.append(
#                         '>' + line[1:6] + '\n' +
#                         next(sequence_file).replace('*',  '').replace('X',  ''))
#                     # Remove the number of remaining IDs by one.
#                     ids_count -= 1
#                     # If all IDs have been found once,  end the search.
#                     if ids_count == 0:
#                         break
#     except IOError:
#         print_error('Cannot find OMA orthologs sequences. '
#                     'The OMA sequence file does not exist!')
#         sys.exit()
#
#     # Fill the orthologous sequences file with the sequences of the
#     # found orthologs. Since newline characters were not trimmed during
#     # the search,  the contents are written down without a delimiter.
#     fill_file(orth_file,  "".join(orth_file_content))


# Helper functions and classes for general use. #
# Splitting lines into columns


def trim_line(line):
    """ Removes the newline at the end of the line. """

    return line.split("\n")[0]


def split_line(line):
    """ Removes the newline at the end of the line and splits the line
    by tabs. """

    return trim_line(line).split("\t")


# File handling

def check_file_exists(filename):
    """ Checks whether the file exists. """
    return os.path.exists(filename)


def fill_file(target,  contents):
    """ Fills a file with the contents. The file is deleted when the
    user interrupts the process. This function mainly prevents partially
     filled zombie files. """

    try:
        with open(target,  'w') as target_file:
            target_file.write(contents)
    except KeyboardInterrupt:
        os.system('rm {0}'.format(target_file))
        sys.exit()

    # removes the variable that contains the target file contents.
    del contents


def generate_lines(target):
    """ Iterates a target file line by line and functions as a generator
    of lines. """

    with open(target,  'r') as target_file:
        for line in target_file:
            yield line


def generate_splitted_lines(target):
    """ Iterates a target file line by line and functions as a generator
    of lines. """

    for line in generate_lines(target):
        yield split_line(line)


# Wrapping a standard routine for every major function call.

def exec_step(
        permission,  output_file,  found_in_cache_message,
        running_function_message,  timestop_message, cache, called_function,
        *args):
    """ Calls a function after checking the config whether this
    function is allowed to run. Then it ensures that no previous result
    is present in cache. The function call is time stopped for passed
    processing time. """

    # Checks the configuration file whether this function is supposed to
    # be called.
    if permission:
        # If a previous result has been cached already, skip the
        # function call and inform the user.
        if (cache
                and output_file is not None
                and cache.in_cache(output_file)):
            print_progress(found_in_cache_message)
        else:
            # Inform the user that the function is called. This is done
            # before stopping the time.
            print_progress(running_function_message)
            # Start stopping the process time. It requires a
            # preexisting instance of TimeStop because the species_id
            # is stored as a class property.
            timer = time_report()
            # Call the actual function. Pass positional arguments.
            called_function(*args)
            # Inform the user about the passed process time.
            timer.print_time(timestop_message)

            # If the function is requested,  but it produces no output,
            # throw an error and exit the program.
            if output_file is not None and not output_file.exists():
                print_error(f'{running_function_message} did not produce an '
                            'output file')
                sys.exit()


class cache_manager:
    """ Saves the reuse_cache and cache_dir config parameter to access
    the cache. """
    __slots__ = ['reuse_cache', 'cache_dir']

    def __init__(self,  prot_config):
        """ Class initialiser """
        self.reuse_cache = prot_config.reuse_cache
        self.cache_dir = prot_config.path_cache

    def in_cache(self, filename, ignore_reuse_cache=False):
        """ Checks if the filename is inside the cache directory.
        If it is,  create a symbolic link inside the current directory."""

        if ignore_reuse_cache:
            use_cache = True
        else:
            use_cache = self.reuse_cache

        cache_path = self.cache_dir / filename.name

        if (use_cache and cache_path.exists()):
            print_progress(f'{filename} found in cache directory. Reusing it.')
            cache_path.symlink_to(Path.cwd())
            return True
        return False

    def copy_to_cache(self,  filename):
        """ Copies the file into the cache directory """
        os.system('cp -a {0} {1}'.format(filename,  self.cache_dir + "/" + filename))
