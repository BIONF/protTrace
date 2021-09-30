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
from subprocess import (
    run as subprocess_run,
    DEVNULL,
    STDOUT,
    PIPE,
    CalledProcessError,
    check_output,
    Popen
)

from dendropy import TreeList
from multiprocessing import Pool
from pathlib import Path

from utils.data_api import proteome
from utils.file import dir_move
from utils.log import (
    print_progress,
    print_warning,
    print_error,
    time_report
)

# Module where actual REvolver / BLAST cycles run takes place
# After every REvolver run, detection probability is checked with reciprocal
# Blast search and stored in a hash table. Using the final hash table, decay
# rate is calculated using decay script


def main(p_id, config):

    # We only need the main ID from the prot_id class
    prot_id = p_id.id

    prot_config = config

    cache = prot_config.reuse_cache
    nr_proc = prot_config.nr_processors

    work_dir = dir_move(prot_config.path_work_dir / prot_id)

    xml_file = Path('revolver_config_' + prot_id + '.xml')

    """ Prepare the query species proteome as a BLAST database within the
    current working directory. """

    # Instantiates the proteome of the query species
    q_proteome = proteome(prot_config, prot_config.species)

    # Creates the filename for the link we use in this module.
    proteome_file = work_dir.current / ('proteome_' + prot_id)

    # The proteome instance asks his source api to link the original proteome
    # file. With a symlink, we need to make sure that the link does not exist
    # yet.
    if not proteome_file.exists():
        q_proteome.transfer_proteome(proteome_file)

    # Runs Makeblastdb
    prepare_blastdb(prot_config, str(proteome_file))

    """ Prepare the simulation tree. """

    trees = TreeList.get_from_path(prot_config.simulation_tree, 'newick')
    taxonset = []
    for element in trees.taxon_namespace:
        taxonset.append(str(element).replace("'", ""))
    taxonset = taxonset[::-1]

    def gen_REvolver_args(prot_id, prot_config, xml_file, proteome_file,
                          taxonset):
        for run in range(prot_config.simulation_runs):
            yield [run, prot_id, prot_config, xml_file, proteome_file,
                   taxonset]

    def run_traceability(prot_id, prot_config, xml_file, proteome_file,
                         taxonset):

        print_progress('Running REvolver / BLAST cycles')
        timestop = time_report()
        results = None
        try:
            with Pool(processes=nr_proc) as pool:
                results = pool.map(actual_traceability_calculation,
                                   gen_REvolver_args(prot_id, prot_config,
                                                     xml_file, proteome_file,
                                                     taxonset))
        except KeyboardInterrupt:
            sys.exit('The user interrupted REvolver.')
        except Exception as e:
            raise e
            print_warning('An error occured during the self identification '
                          'process')
            pass

        if results is None:
            # Make sure to fall back to default values!
            return

        timestop.print_time('REvolver/BLAST', 'mins')

        with open(f'full_decay_results_{prot_id}.txt', 'w') as ffull,\
                open(f'decay_summary_{prot_id}.txt', 'w') as fsum:

            detection_probability = {}
            for res in results:
                for key, value in res.items():
                    if key not in detection_probability.keys():
                        detection_probability[key] = []
                        detection_probability[key].append(value)
                    else:
                        detection_probability[key].append(value)
            for taxa in taxonset:
                ffull.write(taxa + ' ')
                count = 0
                for element in detection_probability[taxa]:
                    ffull.write(str(element))
                    count += int(element)
                ffull.write('\n')
                fsum.write(str(float(count) / float(
                    prot_config.simulation_runs)) + '\n')

        print_progress('Calculating decay parameters.')
        decayParams(prot_config.R, prot_id, prot_config.decay_script)

    if cache and Path(f'decay_summary_{prot_id}.txt_parameter').exists():
        pass
    else:
        run_traceability(prot_id, prot_config, xml_file, proteome_file,
                         taxonset)

    # To remove this temporary file, we need access to the query prot_id class.
    if prot_config.delete_temp:
        p_id.remove_fasta()
    work_dir.reverse()


def prepare_blastdb(prot_config, proteome_file):
    """ Creating a BLAST database for future searches. """
    print_progress('Making a BLAST database of the gene set.')
    cmd = [prot_config.makeblastdb, '-in', proteome_file, '-input_type',
           'fasta', '-dbtype', 'prot']
    try:
        subprocess_run(cmd, check=True, stdout=DEVNULL)
    except CalledProcessError:
        print_error('The query gene set could not be processed into a BLASTDB')

    # Legacy command:
    # os.system('{0} -in {1} -input_type fasta -dbtype prot'
    #           .format(config.path_makeblastdb,  proteome_file))


def actual_traceability_calculation(args):
    """ Perform a single REvolver into BLAST run. """

    run = args[0]
    prot_id = args[1]
    prot_config = args[2]
    xml_file = args[3]
    proteome_file = args[4]
    taxonset = args[5]

    temp_revolver_config_file = Path(f'{run}_revolver_config.xml')
    os.system(f'cp {str(xml_file)} {str(temp_revolver_config_file)}')
    temp_data = temp_revolver_config_file.open('r').read()

    with temp_revolver_config_file.open('w') as ftemp:
        ftemp.write(temp_data.replace('REvolver_output',
                                      'REvolver_output_' + str(run)))

    revolver_output_dir = Path(f'REvolver_output_{run}')
    if not revolver_output_dir.exists():
        revolver_output_dir.mkdir()

    detection_probability = {}

    print_progress(f'Run: {run}')

    success = False
    trials = 0
    while(not success and trials < 10):
        trials += 1
        try:
            try:
                run_revolver(temp_revolver_config_file, prot_config)
            except CalledProcessError as e:
                raise e
                # The exception is raised by the subprocess routine that
                # calles REvolver
                print_warning('REvolver threw an exception!')
                break
            try:
                blastOutput = run_blast(prot_config.blastp, prot_id,
                                        proteome_file, revolver_output_dir)
            except CalledProcessError as e:
                raise e
                # The exception is raised by the subprocess routine that
                # calles BLASTP
                print_warning('BLASTP threw an exception during the reblast!')
                break
            for taxa in taxonset:
                detection = 0
                linecount = 1
                for line in blastOutput.split('\n'):
                    if taxa == line.split('\t')[0]:
                        if line.split('\t')[1] == prot_id:
                            # blastHitId:
                            detection = 1
                            break
                        linecount += 1
                        if linecount > 5:
                            break
                detection_probability[taxa] = detection
            success = True
        except KeyboardInterrupt:
            print('Keyboard interruption by user!')
            raise Exception

    if trials >= 10:
        print_warning('Too many trials for REvolver!')

    if prot_config.delete_temp:
        temp_revolver_config_file.unlink()
        os.system(f'rm -r {revolver_output_dir}')

    return detection_probability


def decayParams(r, prot_id, decay_script):
    """ Runs the Rscript to calculate the decay parameters from the REvolver
    runs. """
    # command = (f'{r} --quiet --vanilla {decay_script} '
    #            f'decay_summary_{prot_id}.txt')

    output_file = f'decay_summary_{prot_id}.txt'
    command = [r, '--quiet', '--vanilla', decay_script, output_file]

    # if __debug__:
    #     print(f'Decay parameter calculation command: {command}')

    subprocess_run(command, check=True)


def run_revolver(xml_file, prot_config):
    command = [prot_config.java, '-Xmx2G', '-Xms2G', '-cp',
               prot_config.REvolver, 'revolver', str(xml_file)]
    try:
        check_output(command, stderr=STDOUT)
    except CalledProcessError as e:
        # The exception gets caught to not halt the entire ProtTrace process.
        # This is especially important when processing multiple query proteins.
        # Since each protein is processed consecutively.

        # Any exception, that originates from deleting the entire sequence
        # and trying to evolving it further with feature constraints,
        # is reduced to few sentences to not clutter the log.
        # For this, the message is splitted for each sentence and then stripped
        # off the sometimes varying line numbers in the java sourcecode.
        message_lines = set([m.split("(")[0] for m
                             in e.output.decode('utf-8').split("\n")])
        if '	at controller.seqGeneration.HmmSeqEvolver.'\
                'evolveSequence_waitingTime' in message_lines:
            print_warning('The feature-dependent HmmSeqEvolver cannot '
                          'simulate the evolution/deletion/insertion waiting '
                          'time!\nThe sequence evolved most likely at rapid '
                          'speeds and got fully deleted!')
            pass
        elif ('For input string: \">>\"') in message_lines:
            print_warning('REvolver had an issue with a wrongly formatted '
                          f'input sequence in run {xml_file}')
            pass
        else:
            print(e.output.decode('utf-8'))
            # This function needs to raise the exception, otherwise the
            # wrapping functions try-catch will not detect the exception!
            pass
            # raise e


def run_blast(blastp, prot_id, proteome, revolverOut):
    # Performs the reblast against the query proteome.

    # Declares the BLASTP process.
    blast_process = Popen([blastp, '-query', f'{revolverOut}/out.fa', '-db',
                           proteome, '-outfmt', '6 qseqid sseqid'],
                          stdout=PIPE, stderr=PIPE)
    # Waits for BLASTP to be finished and retrieves the result and error
    # messages.
    result, error = blast_process.communicate()
    # Necessary, since in python 3, subprocess outputs bytes. In python 3,
    # bytes and strings are not interchangeable anymore.
    result = result.decode('utf-8')
    error = error.decode('utf-8')
    # This error message is caught and printed only once.
    if error[0:54] == 'BLAST engine error: Warning: Sequence contains no data':
        print('BLAST engine error: Warning: Sequence contains no data')
        print('Either an sequence in the simulation step was deleted to zero '
              'amino-acids or REvolver did not execute properly!')
    elif error == 'Command line argument error: Argument \"query\". File is '\
            'not accessible: {0}/out.fa\n'.format(revolverOut):
        print('REvolver did not produce an output fasta file!')
    elif error != "":
        print(error)
    return result
