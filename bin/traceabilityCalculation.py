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
import configure
import subprocess
import dendropy
import time
from multiprocessing import Pool

### Module where actual REvolver / BLAST cycles run takes place
### After every REvolver run, detection probability is checked with reciprocal Blast search and stored in a hash table
### Using the final hash table, decay rate is calculated using decay script

def main(p_id, config_file):
    # Global variables are needed, since the multiprocessing.Pool.map function only maps one parameter onto each run.
    # All other outside parameters in use need to be globally available.
    global prot_id
    prot_id = p_id

    global prot_config
    prot_config = configure.setParams(config_file)
    cache = prot_config.reuse_cache
    nr_proc = prot_config.nr_processors

    rootDir = os.getcwd()
    work_dir = prot_config.path_work_dir + '/' + prot_id
    os.chdir(work_dir)

    if os.path.exists(work_dir + '/omaId.txt'):
        global blastHitId
        blastHitId = open(work_dir + '/omaId.txt').read().split('\n')[0]
    else:
        sys.exit('### ERROR: No reciprocal BLAST hit id found!!!! ###')

    global xml_file
    xml_file = 'revolver_config_' + prot_id + '.xml'

    global proteome_file
    proteome_file = work_dir + '/' + 'proteome_' + prot_id

    trees = dendropy.TreeList.get_from_path(prot_config.simulation_tree, "newick")
    global taxonset
    taxonset = []
    for element in trees.taxon_namespace:
        taxonset.append(str(element).replace("'", ""))
    taxonset = taxonset[::-1]

    print('##### Running REvolver / BLAST cycles: #####')
    start_time = time.time()

    if cache and os.path.exists('decay_summary_%s.txt_parameter' %prot_id):
        pass
    else:
        try:
            pool = Pool(processes=nr_proc)
            results = pool.map(actual_traceability_calculation, range(prot_config.simulation_runs))
        except KeyboardInterrupt:
            pool.terminate()
            pool.close()
            sys.exit("Interrupting REvolver")
        #except:
        #	print("ERROR: Multiprocessing step <-> Traceability Calculations.")
        #	pass
        finally:
            pool.terminate()
            pool.join()

        print('#####\tTIME TAKEN: %s mins REvolver/BLAST#####' %((time.time() - start_time) / 60))

        ffull = open('full_decay_results_%s.txt' %prot_id, 'w')
        fsum = open('decay_summary_%s.txt' %prot_id, 'w')

        detection_probability = {}
        for res in results:
            for key, value in res.items():
                if not key in detection_probability.keys():
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
            fsum.write(str(float(count) / float(prot_config.simulation_runs)) + '\n')
        ffull.close()
        fsum.close()

        print('##### Calculating decay parameters #####')
        decayParams(prot_config.R, prot_id, prot_config.decay_script)

    os.chdir(rootDir)

def actual_traceability_calculation(run):
    temp_revolver_config_file = "%s_revolver_config.xml" %run
    os.system('cp %s %s' %(xml_file,temp_revolver_config_file))
    temp_data = open(temp_revolver_config_file).read()
    ftemp = open(temp_revolver_config_file, 'w')
    ftemp.write(temp_data.replace('REvolver_output', 'REvolver_output_' + str(run)))
    ftemp.close()

    revolver_output_dir = "REvolver_output_%s" %run
    if not os.path.exists(revolver_output_dir):
        os.mkdir(revolver_output_dir)

    detection_probability = {}

    print('Run: {0}'.format(run))

    success = False
    trials = 0
    while(not success and trials < 10):
        trials += 1
        try:
            try:
                run_revolver(prot_config.REvolver, temp_revolver_config_file)
            except Exception as e:
                # The exception is raised by the subprocess routine that calles REvolver
                print("REvolver threw an exception!")
                break
            try:
                blastOutput = run_blast(prot_config.blastp, prot_id, proteome_file, revolver_output_dir)
            except Exception as e:
                # The exception is raised by the subprocess routine that calles BLASTP
                print("BLASTP threw an exception during the reblast!")
                break
            for taxa in taxonset:
                detection = 0
                linecount = 1
                for line in blastOutput.split('\n'):
                    #print(line)
                    if taxa == line.split('\t')[0]:
                        if line.split('\t')[1] == blastHitId:
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
        print('TOO MANY TRIALS FOR REVOLVER!!! Check REvolver configuration file.')

    if prot_config.delete_temp:
        os.system('rm -rf %s' %revolver_output_dir)
        os.system('rm -rf %s' %temp_revolver_config_file)

    return detection_probability

def decayParams(r, prot_id, decay_script):
    command = '%s --quiet --vanilla %s decay_summary_%s.txt' %(r, decay_script, prot_id)
    print('##### Decay parameter calculation command: ', command)
    os.system(command)

def run_revolver(REvolver, xml_file):
    command = ["java","-Xmx2G","-Xms2G","-cp",REvolver,"revolver",xml_file]
    try:
        subprocess.check_output(command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        # The exception gets caught to not halt the entire ProtTrace process.
        # This is especially important when processing multiple query proteins.
        # Since each protein is processed consecutively.

        # Any exception, that originates from deleting the entire sequence
        # and trying to evolving it further with feature constraints,
        # is reduced to few sentences to not clutter the log.
        # For this, the message is splitted for each sentence and then stripped off the
        # sometimes varying line numbers in the java sourcecode.
        message_lines = set([m.split("(")[0] for m in e.output.decode('utf-8').split("\n")])
        if "	at controller.seqGeneration.HmmSeqEvolver.evolveSequence_waitingTime" in message_lines:
            print("The feature-dependent HmmSeqEvolver cannot simulate the evolution/deletion/insertion waiting time!\n\
The sequence evolved most likely at rapid speeds and got fully deleted!")
            pass
        else:
            print(e.output.decode('utf-8'))
            raise e

def run_blast(blastp, prot_id, proteome, revolverOut):
    # Performs the reblast against the query proteome.

    # Declares the BLASTP process.
    blast_process = subprocess.Popen([blastp,"-query","{0}/out.fa".format(revolverOut),"-db",proteome,"-outfmt","6 qseqid sseqid"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # Waits for BLASTP to be finished and retrieves the result and error messages.
    result, error = blast_process.communicate()
    # Necessary, since in python 3, subprocess outputs bytes. In python 3, bytes and strings are not interchangeable anymore.
    result = result.decode('utf-8')
    error = error.decode('utf-8')
    # This error message is caught and print(ed only once.)
    if error[0:54] == "BLAST engine error: Warning: Sequence contains no data":
        print("BLAST engine error: Warning: Sequence contains no data")
        print("Either an sequence in the simulation step was deleted to zero amino-acids or REvolver did not execute properly!")
    elif error == "Command line argument error: Argument \"query\". File is not accessible:  `{0}/out.fa'\n".format(revolverOut):
        print("REvolver did not produce an output fasta file!")
    elif error != "":
        print(error)
    return result

