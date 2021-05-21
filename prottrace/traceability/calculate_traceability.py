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

import math
from pathlib import Path

from Bio import SeqIO

from utils.data_api import (
    species_mapping,
    get_species_distance
)
from utils.file import dir_move
from utils.log import print_progress

# Script to colourize the species according to the traceability values
# Reads in the generated nexus file by PTP and edits it according to our need
# Represents the traceabilities as different colours on the tip labels
# INPUT:
# 1. Nexus tree file
# 2. HaMStR mapping file
# 3. Species Id
# 4. Species tree
# 5. Decay result file


# Read the pairwise species maximum likelihood distance from
# the table file or from cache.
def gather_max_lik_dist(species1, species2, distance_matrix_file, cache):

    # To check if the species are present in the species ML dist
    # matrix file. Otherwise, parse likelihood from cache directory.
    not_found_ml1 = True
    not_found_ml2 = True

    # Read in the distance matrix
    species_ml = distance_matrix_file.open('r').read().split('\n')

    for line in range(len(species_ml) - 1):
        if species_ml[line].split('\t')[0] == species1:
            rowIndex = line
            not_found_ml1 = False
        elif species_ml[line].split('\t')[0] == species2:
            columnIndex = line
            not_found_ml2 = False

    ml_1 = cache / f'{species1}_{species2}.lik'
    ml_2 = cache / f'{species2}_{species1}.lik'

    if not_found_ml1 or not_found_ml2:
        # Checking for the likelihood score in cache directory´
        if ml_1.exists():
            return float(ml_1.open('r').read().split('\n')[0])
        elif ml_2.exists():
            return float(ml_2.open('r').read().split('\n')[0])
        else:
            return 1.00
    else:
        if not species_ml[rowIndex].split('\t')[columnIndex] == "NA":
            return float(species_ml[rowIndex].split('\t')[columnIndex])
        else:
            return 1.00


# def calculate_traceability(query, subject, decay_rate, decay_pop,
#                            species_dist_file, cache):
def calculate_traceability(query, subject, decay_rate, decay_pop, config):

    # Do not bother to collect the ML distance for trivial cases.
    if decay_rate < 0.01 or query.species_id == subject.species_id:
        traceability = 1
    else:
        ml_dist = get_species_distance(config, query, subject, default=1.0)
        # ml_dist = gather_max_lik_dist(query, subject, species_dist_file, cache)

        # This is the calculation of the protein's traceability index in a
        # species (specified by ml_dist).
        traceability = 1 - ((decay_pop * math.exp(decay_rate * ml_dist)) /
                            (1 + decay_pop *
                             (math.exp(decay_rate * ml_dist) - 1)))

    return traceability


# Writes the traceability values into two table files
# The table can be written into two styles:
# table        - A regular table with target species names
# phyloprofile - A table readable by PhyloProfile that uses NCBI Ids
def create_traceability_output_file(query, mappings, compute_fas,
                                    config):

    # The protein ID
    qu_prot = query.id

    # Retrieve the presence of orthologs among all other species
    orth_species_set = [record.description for record in
                        SeqIO.parse(f'ogSeqs_{qu_prot}.fa', 'fasta')]
    # if orth_file_name.exists():
    #     with orth_file_name.open('r') as of:
    #         # Distills the sequence headers to a set of species with orthologs
    #         # If the line contains underscores to differentiate same-species
    #         # sequences, assume the part before the first underscore to be the
    #         # species id. The query species is also listed in the orthologs
    #         # file.
    #         orth_species_set = {line[1:].split('_')[1] if '_' in line
    #                             else line[1:] for line in of if line[0] == '>'}

    # # FAS scores are only retrieved if the option was turned on in the config
    # if compute_fas:
    #     fas_file_name = Path(f'ogSeqs_{qu_prot}.fasScore')
    #     if fas_file_name.exists():
    #         with fas_file_name.open('r') as ff:
    #             # Creates a dictionary of IDs to FAS scores
    #             # The ID resolution is the long beginning part of the line
    #             fas_scores = {line.split()[2].split('_')[1] if '_' in line
    #                           else line.split()[2]: line.split()[-1]
    #                           for line in ff}

    # We put all output table lines into this list
    output_lines = []
    # We put all output phyloprofile-compatible table lines into this list
    output_phyloprofile_lines = []

    # The phyloprofile style includes an initial line with column headers
    # if compute_fas:
    #     output_phyloprofile_lines.append(["geneID", "ncbiID", "orthID",
    #                                       "FAS", "Traceability"])
    # else:
    output_phyloprofile_lines.append(["geneID", "ncbiID", "orthID",
                                      "Traceability"])

    # Get the decay rate and decay pop of the protein to calculate the
    # traceability
    with Path(f'decay_summary_{qu_prot}.txt_parameter').open('r') as t_decay_f:
        written_parameters = t_decay_f.read().split("\n")
        decay_rate = float(written_parameters[0].rstrip())
        decay_pop = float(written_parameters[1].rstrip())

    # Fill the output table lines with information.
    # Start with the species names from the mapping file.

    for spec in mappings.mappings:
        trace = calculate_traceability(query.spec_names, spec, decay_rate,
                                       decay_pop, config)
        # Add lines to the table output
        output_lines.append([query.spec_name, spec.taxon, str(trace)])

        # Add lines to the phyloprofile table output. The presence and absence
        # of orthologs is expressed in binary.
        orthology = 0
        if spec.species_id in orth_species_set:
            orthology = 1
        # FAS computation is no longer supported in ProtTrace. Choose
        # dedicated software instead.
        # if compute_fas:
        #     fas = 'NA'
        #     if species_id == qu_prot:
        #         fas = 1
        #     if species_id != qu_prot:
        #         fas = fas_scores[species_id]
        #     output_phyloprofile_lines.append([qu_prot,
        #                                       "ncbi" + species_line[2],
        #                                       str(orthology),
        #                                       str(fas),
        #                                       trace])
        # else:
        output_phyloprofile_lines.append([qu_prot, f'ncbi{spec.ncbi}',
                                          str(orthology), str(trace)])

    # Compile the output lines into one string with newlines
    # and the intended column separator
    output = '\n'.join(['\t'.join(columns) for columns in output_lines])
    output_phyloprofile = '\n'.join(['\t'.join(columns) for columns in
                                     output_phyloprofile_lines])

    print_progress('Writing the traceabilites into output files.')

    # Write both output files
    with Path(f'trace_results_{qu_prot}.txt').open('w') as output_file:
        output_file.write(output)
    with Path(f'{qu_prot}_phyloMatrix.txt').open('w') as output_pp_file:
        output_pp_file.write(output_phyloprofile)


# def main(query_id, cache, nexusTreeFile, mapFile, protId, spTree,
# plotFigTree,
#          speciesMaxLikFile, speciesId, cache_dir, intended_fas_score_calc,
# config):
def main(query, config):

    # Assuming that we did not need the species mapping for a long time,
    # we can create it again from the species mapping file, rather than
    # occupying the extra RAM all the time.
    spec_mapping = species_mapping(config)
    # cache = config.path_cache
    compute_fas = config.fas_score
    # species_dist_file = config.ML_table
    matrixDict = {}
    fas_file = Path(f'ogSeqs_{query.id}.fasScore')
    orth_file = Path(f'ogSeqs_{query.id}.fa')
    # decayRate = 0.1
    # decayPop = 0.04

    work_dir = dir_move(config.path_work_dir / query.id)

    print_progress('Creating the traceability table.')

    # Iterate through all species in the species mapping file
    for unit in spec_mapping.all_species():
        orth = "NA"
        fas = "NA"
        if unit not in matrixDict:
            matrixDict[unit] = []

        # The query species gets FAS score and orthology presence value of 1
        if unit == query.id:
            if fas_file.exists():
                with fas_file.open('r') as fas_f:
                    fas_entries = fas_f.read().split('\n')[0]
                    if fas_entries != "":
                        orth = fas_entries.split()[1]
                    fas = "1.00"
            else:
                orth = "1"

        # Any species that is not the query species
        else:
            # Adds FAS score information to the output table
            if fas_file.exists():
                with fas_file.open('r') as fas_f:
                    for line2 in fas_f:
                        line2_species_id = (line2.split()[2].split('_')[1]
                                            if '_' in line2
                                            else line2.split()[2])
                        if line2_species_id == unit:
                            orth = line2.split()[3]
                            fas = line2.split()[-1]
            # Adds orthology information to the output table
            elif orth_file.exists():
                orth = "0"
                with orth_file.open('r') as orth_f:
                    for line2 in orth_f:
                        if '>' in line2:
                            line2_species_id = (line2.split()[0].split('_')[1]
                                                if '_' in line2
                                                else line2.split()[0][1:])
                            if line2_species_id == unit:
                                orth = "1"

        # Whatever orth and fas values were set above, now they are added for
        # this species.
        matrixDict[unit].append(orth + '#' + fas)

    # print('Dictionary length: ', len(matrixDict))
    # try:
    #     fdogMapFile = open(mapFile).read().split('\n')
    #     speciesTree = open(spTree).read()
#   #      decayRate = float(open('decay_summary_%s.txt_parameter' %protId)
# .r# ead().split('\n')[1])
#   #      decayPop = float(open('decay_summary_%s.txt_parameter' %protId)
# .r# ead().split('\n')[0])

    # except IOError:
    #     print('ERROR: Colourizing tree encountered problem!!!')

    # Calculates and writes the protein traceability into output files
    create_traceability_output_file(query, spec_mapping, compute_fas, config)

    # Write the output table file
#    traceResults = open('trace_results_%s.txt' %protId, 'w')

    work_dir.reverse()

#    if os.path.exists(nexusTreeFile):
#    #speciesId = nexusTreeFile.split('_')[1].split('.')[0][:5]
#        for i in range(len(fdogMapFile) - 1):
#            if speciesId == fdogMapFile[i].split('\t')[-1]:
#                speciesName = fdogMapFile[i].split('\t')[1]
#                break
#
#        colourize(speciesName, speciesId, nexusTreeFile, speciesTree, decayRate, decayPop, traceResults, matrixDict)
#
#    traceResults.close()
#
#        # Visualizes the species tree that has been colourized by the
#        # traceability index of the query protein. It outputs the tree
#        # in PDF format.
#        try:
#            # LEGACY code
#            # subprocess.check_output('java -cp %s figtreepdf %s' %(plotFigTree, nexusTreeFile.replace('.nexus', '_edit.nexus')),shell=True)
#            subprocess.check_output(['java','-Djava.awt.headless=true','-cp',plotFigTree,'figtreepdf',nexusTreeFile.replace('.nexus', '_edit.nexus')])
#            print("Done!")
#        except subprocess.CalledProcessError as e:
#            print(e.output)
#            print('WARNING: No representation of traceabilities on tree possible.\nJAVA program figtreepdf not responding!!!')
#
#    # Writing the traceability output as a matrix which can be read with PhyloProfile
#    print('Creating matrix file for PhyloProfile...')
#    matrixFile = open('%s_phyloMatrix.txt' %protId, 'w')
#    if intended_fas_score_calc:
#        matrixFile.write('geneID\tncbiID\torthoID\tFAS\tTraceability\n')
#    else:
#        matrixFile.write('geneID\tncbiID\torthoID\tTraceability\n')
#    for j in range(len(fdogMapFile) - 1):
#        ncbiId = "ncbi" + fdogMapFile[j].split()[-2]
#        oma_id = fdogMapFile[j].split()[-1]
#        for elements in matrixDict[oma_id]:
#            if intended_fas_score_calc:
#                matrixFile.write(protId + '\t' + ncbiId + '\t' + elements.split('#')[0] + '\t' + elements.split('#')[1] + '\t' + elements.split('#')[2] + '\n')
#            else:
#                matrixFile.write(protId + '\t' + ncbiId + '\t' + elements.split('#')[0] + '\t' + elements.split('#')[2] + '\n')
#    matrixFile.close()
