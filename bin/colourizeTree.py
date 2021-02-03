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
import math
import subprocess

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
def gather_max_lik_dist(species1, species2, distanceMatrixFile):

    # To check if the species are present in the species ML dist
    # matrix file. Otherwise, parse likelihood from cache directory.
    flag1 = True
    flag2 = True

    # Read in the distance matrix
    speciesMaxFile = open(distanceMatrixFile).read().split('\n')

    for line in range(len(speciesMaxFile) - 1):
        if speciesMaxFile[line].split('\t')[0] == species1:
            rowIndex = line
            flag1 = False
        elif speciesMaxFile[line].split('\t')[0] == species2:
            columnIndex = line
            flag2 = False

    if flag1 or flag2:
#        print('Checkpoint 2 crossed')
        # Checking for the likelihood score in cache directory
        if os.path.exists(cacheDir + '/' + species1 + '_' + species2 + '.lik'):
            return float(open(cacheDir + '/' + species1 + '_' + species2 + '.lik').read().split('\n')[0])
        elif os.path.exists(cacheDir + '/' + species2 + '_' + species1 + '.lik'):
            return float(open(cacheDir + '/' + species2 + '_' + species1 + '.lik').read().split('\n')[0])
        else:
#            print('No likelihood distance found between species: %s and %s. Using default likelihood distance of 1.0!' %(species1, species2))
            return 1.00
    else:
        if not speciesMaxFile[rowIndex].split('\t')[columnIndex] == "NA":
            return float(speciesMaxFile[rowIndex].split('\t')[columnIndex])
        else:
#            print('No likelihood distance found between species: %s and %s. Using default likelihood distance of 1.0!' %(species1, species2))
            return 1.00

def colourize(speciesName, speciesId, nexusTreeFile, speciesTree, decayRate, decayPop, traceResults, matrixDict):
    tree = open(nexusTreeFile).read().split('\n')
    for i in range(len(tree) - 1):
        if tree[i] == "\ttaxlabels":
            startTaxa = i + 1
        if tree[i] == ";":
            stopTaxa = i
        if tree[i] == "begin trees;":
            treeLine = i + 2
            break
    fnew = open(nexusTreeFile.replace('.nexus', '_edit.nexus'), 'w')
    for i in range(startTaxa):
        fnew.write(tree[i] + '\n')

    for i in range(startTaxa, stopTaxa):
        tempSpecies = tree[i].split()[-1]
        if not tempSpecies == speciesName:
            try:
                colourCode, traceValue = getColourCode(speciesName, tempSpecies, decayRate, decayPop)
                omaName = "NA"

                for j in range(len(hamstrMapFile) - 1):
                    if tempSpecies == hamstrMapFile[j].split('\t')[1]:
                        omaName = hamstrMapFile[j].split('\t')[-1]
                        break
                for elements in matrixDict[omaName]:
                    newElement = elements + '#' + str(traceValue)
                    matrixDict[omaName].remove(elements)
                    matrixDict[omaName].insert(0, newElement)
                traceResults.write(speciesName + '\t' + tempSpecies + '\t' + str(traceValue) +'\n')
                fnew.write(tree[i] + '[&!color=#-' + colourCode + ']' + '\n')
            except:
                sys.exit("ERROR: Check species %s in species tree and used mapping files" %tempSpecies)
                pass

        else:
            fnew.write(tree[i] + '\n')
            for elements in matrixDict[speciesId]:
                newElement = elements + "#1"
                matrixDict[speciesId].remove(elements)
                matrixDict[speciesId].insert(0,newElement)

    for i in range(stopTaxa, treeLine):
        fnew.write(tree[i] + '\n')
    fnew.write(speciesTree.replace('\n', '') + '\n')
    for i in range(treeLine + 1, len(tree)):
        if tree[i] == '	set branchLabels.isShown=true;':
            fnew.write("	set branchLabels.isShown=false;" + '\n')
        else:
            fnew.write(tree[i] + '\n')
    fnew.close()

def getColourCode(spName, tempName, decayRate, decayPop, species_distances_file):
    """ Colour codes the branches of the visualization tree by the
    traceability index. """
    mlDist = gather_max_lik_dist(spName, tempName, species_distances_file)
    #print(mlDist)
    if decayRate < 0.01:
        traceability = 1
    else:
        # This is the calculation of the protein's traceability index in a species (specified by mlDist).
        traceability = 1 - ((decayPop * math.exp(decayRate * mlDist)) / (1 + decayPop * (math.exp(decayRate * mlDist) - 1)))
    #print(traceability)
    if traceability <= 1 and traceability >= 0.9:
        colCode = '16711936'
    elif traceability < 0.9 and traceability >= 0.8:
        colCode = '3604736'
    elif traceability < 0.8 and traceability >= 0.7:
        colCode = '2294016'
    elif traceability < 0.7 and traceability >= 0.6:
        colCode = '983296'
    elif traceability < 0.6 and traceability >= 0.5:
        colCode = '256'
    elif traceability < 0.5:
        colCode = '65536'

    return colCode, traceability


def calculate_traceability(query_species, species_id, decay_rate, decay_pop,
                           species_distances_file):

    ml_dist = gather_max_lik_dist(query_species, species_id,
                                  species_distances_file)
    if decay_rate < 0.01:
        traceability = 1
    else:
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
def create_traceability_output_file(query_species, query_protein,
                                    species_mapping_file, compute_fas,
                                    species_distances_file):

    # Retrieve the presence of orthologs among all other species
    orth_file_name = 'ogSeqs_' + query_protein + '.fa'
    if os.path.exists(orth_file_name):
        with open(orth_file_name) as of:
            # Distills the sequence headers to a set of species with orthologs
            # If the line contains underscores to differentiate same-species
            # sequences, assume the part before the first underscore to be the
            # species id. The query species is also listed in the orthologs
            # file.
            orth_species_set = {line[1:].split('_')[1] if '_' in line
                                else line[1:] for line in of if line[0] == '>'}

    # FAS scores are only retrieved if the option was turned on in the config
    if compute_fas:
        fas_file_name = 'ogSeqs_' + query_protein + '.fasScore'
        if os.path.exists(fas_file_name):
            with open(fas_file_name, 'r') as ff:
                # Creates a dictionary of IDs to FAS scores
                # The ID resolution is the long beginning part of the line
                fas_scores = {line.split()[2].split('_')[1] if '_' in line
                              else line.split()[2]: line.split()[-1]
                              for line in ff}

    # We put all output table lines into this list
    output_lines = []
    # We put all output phyloprofile-compatible table lines into this list
    output_phyloprofile_lines = []

    # The phyloprofile style includes an initial line with column headers
    if compute_fas:
        output_phyloprofile_lines.append(["geneID", "ncbiID", "orthID", "FAS",
                                          "Traceability"])
    else:
        output_phyloprofile_lines.append(["geneID", "ncbiID", "orthID",
                                          "Traceability"])

    # The default output table uses the proper species name, which we will
    # find in the mapping file later phyloprofile uses the query protein name,
    # which is given in the parameter already
    query_name = query_protein

    # Get the decay rate and decay pop of the protein to calculate the
    # traceability
    with open('decay_summary_{0}.txt_parameter'
              .format(query_protein)) as traceability_decay_parameters_file:
        written_parameters = (traceability_decay_parameters_file
                              .read().split("\n"))
        decay_rate = float(written_parameters[0])
        decay_pop = float(written_parameters[1])

    # Fill the output table lines with information
    # Start with the species names from the mapping file
    with open(species_mapping_file, 'r') as mapping_file:
        counter = 0
        query_found = 0
        for line in mapping_file:
            # Make line elements retrievable by index
            species_line = line.rstrip().split("\t")

            species_id = species_line[-1]

            # Get the evolutionary traceability in this species before
            # deciding for a style
            species_traceability = calculate_traceability(query_species,
                                        species_id, decay_rate, decay_pop,
                                        species_distances_file)

            # Add lines to the table output
            output_lines.append([query_name, species_line[1],
                                 str(species_traceability)])

            # Add lines to the phyloprofile table output
            orthology = 0
            if species_id in orth_species_set:
                orthology = 1
            if compute_fas:
                fas = 'NA'
                if species_id == query_protein:
                    fas = 1
                if species_id != query_protein:
                    fas = fas_scores[species_id]
                output_phyloprofile_lines.append([query_protein,
                                                  "ncbi" + species_line[2],
                                                  str(orthology),
                                                  str(fas),
                                                  species_traceability])
            else:
                output_phyloprofile_lines.append([query_protein,
                                                  "ncbi" + species_line[2],
                                                  str(orthology),
                                                  str(species_traceability)])

            # The counter allows us to modify the query_name later when we
            # find the query's line
            if query_found == 0 and species_line[-1] == query_species:
                query_name = species_line[1]
                # We freeze the counter and replace the first column in the
                # default table output after the loop
                query_found = 1
            else:
                counter += 1

    # The query name has been likely found in the middle of the mapping file
    # Here, we replace the first column of any species we wrote
    # If this step fails, the first column is populated with the query protein
    # ID
    for c in range(counter):
        output_lines[c][0] = query_name

    # Compile the output lines into one string with newlines
    # and the intended column separator
    output = "\n".join(["\t".join(line) for line in output_lines])
    output_phyloprofile = "\n".join(["\t".join(line) for line in
                                     output_phyloprofile_lines])

    # Write both output files
    with open('trace_results_{0}.txt'.format(query_protein),
              'w') as output_file:
        output_file.write(output)
    with open('{0}_phyloMatrix.txt'.format(query_protein),
              'w') as output_pp_file:
        output_pp_file.write(output_phyloprofile)


def main(nexusTreeFile, mapFile, protId, spTree, plotFigTree,
         speciesMaxLikFile, speciesId, cache_dir, intended_fas_score_calc):
    global sp_max_file, hamstrMapFile, cacheDir

    matrixDict = {}
    fasFile = 'ogSeqs_' + protId + '.fasScore'
    orthFile = 'ogSeqs_' + protId + '.fa'
    sp_max_file = speciesMaxLikFile
    decayRate = 0.1
    decayPop = 0.04
    cacheDir = cache_dir

    # Iterate through all species in the species mapping file
    for line in open(mapFile):
        orth = "NA"
        fas = "NA"
        currentSpecies = line.split()[-1]
        if currentSpecies not in matrixDict:
            matrixDict[currentSpecies] = []

        # The query species gets FAS score and orthology presence value of 1
        if currentSpecies == speciesId:
            if os.path.exists(fasFile):
                if open(fasFile).read().split('\n')[0] != "":
                    orth = open(fasFile).read().split('\n')[0].split()[1]
                fas = "1.00"
            else:
                orth = "1"

            matrixDict[currentSpecies].append(orth + '#' + fas)
        # Any species that is not the query species
        else:
            foundSpeciesFlag = False
            # Adds FAS score information to the output table
            if os.path.exists(fasFile):
                for line2 in open(fasFile):
                    line2_species_id = (line2.split()[2].split('_')[1]
                                        if '_' in line2 else line2.split()[2])
                    if line2_species_id == currentSpecies:
                        foundSpeciesFlag = True
                        orth = line2.split()[3]
                        fas = line2.split()[-1]
                        matrixDict[currentSpecies].append(orth + '#' + fas)
            # Adds orthology information to the output table
            elif os.path.exists(orthFile):
                orth = "0"
                for line2 in open(orthFile):
                    if '>' in line2:
                        line2_species_id = (line2.split()[0].split('_')[1]
                                            if '_' in line2
                                            else line2.split()[0][1:])
                        if line2_species_id == currentSpecies:
                            foundSpeciesFlag = True
                            orth = "1"
                            matrixDict[currentSpecies].append(orth + '#' + fas)
            if not foundSpeciesFlag:
                matrixDict[currentSpecies].append(orth + '#' + fas)

    # print('Dictionary length: ', len(matrixDict))
    try:
        hamstrMapFile = open(mapFile).read().split('\n')
        speciesTree = open(spTree).read()
#        decayRate = float(open('decay_summary_%s.txt_parameter' %protId)
# .read().split('\n')[1])
#        decayPop = float(open('decay_summary_%s.txt_parameter' %protId)
# .read().split('\n')[0])

    except IOError:
        print('ERROR: Colourizing tree encountered problem!!!')

    # Calculates and writes the protein traceability into output files
    create_traceability_output_file(speciesId, protId, mapFile,
                                    intended_fas_score_calc, speciesMaxLikFile)

    # Write the output table file
#    traceResults = open('trace_results_%s.txt' %protId, 'w')

#    if os.path.exists(nexusTreeFile):
#    #speciesId = nexusTreeFile.split('_')[1].split('.')[0][:5]
#        for i in range(len(hamstrMapFile) - 1):
#            if speciesId == hamstrMapFile[i].split('\t')[-1]:
#                speciesName = hamstrMapFile[i].split('\t')[1]
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
#    for j in range(len(hamstrMapFile) - 1):
#        ncbiId = "ncbi" + hamstrMapFile[j].split()[-2]
#        oma_id = hamstrMapFile[j].split()[-1]
#        for elements in matrixDict[oma_id]:
#            if intended_fas_score_calc:
#                matrixFile.write(protId + '\t' + ncbiId + '\t' + elements.split('#')[0] + '\t' + elements.split('#')[1] + '\t' + elements.split('#')[2] + '\n')
#            else:
#                matrixFile.write(protId + '\t' + ncbiId + '\t' + elements.split('#')[0] + '\t' + elements.split('#')[2] + '\n')
#    matrixFile.close()

