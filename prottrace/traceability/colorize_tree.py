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

from dendropy import TreeList

from utils.configure import set_params
from utils.log import print_progress
from utils.file import dir_move

""" This module colorizes a newick species tree according to the traceability
    of the query species. """


def main(query, config):

    prot_id = query.id
    species_id = query.spec_id

    work_dir = dir_move(config.path_work_dir / prot_id)

    nexus_file = 'nexus_' + prot_id + '.nexus'

    # Save the species name into a variable taxonset
    trees = TreeList.get_from_path(config.species_tree, "newick")
    taxonset = []
    for element in trees.taxon_namespace:
        taxonset.append(str(element).replace("'", "").replace(" ", "_"))

    generateNexusFile(nexus_file, taxonset, config)
    # colourizeTree.main(nexus_file, config.fdog_oma_tree_map,
    #                    prot_id, config.species_tree,
    #                    config.plot_figtree,
    #                    config.ML_matrix,
    #                    species_id, cache_dir, config.fas_score)

    work_dir.reverse()


def colorize(query_species, speciesId, nexusTreeFile, speciesTree, decayRate,
             decayPop, traceResults, matrixDict):
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

    def gen_species(tree, start, stop):
        """ Iterates through the nexus file and returns ever species name. """
        for i in range(start, stop):
            yield tree[i].split()[-1]

    for taxon in gen_species(tree, startTaxa, stopTaxa):
        if not taxon == query_species:
            try:
                colourCode, traceValue = getColourCode(query_species,
                                                       taxon,
                                                       decayRate, decayPop)
                omaName = "NA"

                for j in range(len(fdogMapFile) - 1):
                    if taxon == fdogMapFile[j].split('\t')[1]:
                        omaName = fdogMapFile[j].split('\t')[-1]
                        break
                for elements in matrixDict[omaName]:
                    newElement = elements + '#' + str(traceValue)
                    matrixDict[omaName].remove(elements)
                    matrixDict[omaName].insert(0, newElement)
                traceResults.write(species_name + '\t' + taxon + '\t' +
                                   str(traceValue) +'\n')
                fnew.write(tree[i] + '[&!color=#-' + colourCode + ']' + '\n')
            except:
                print_error(f'Check species {taxon} in species tree '
                            'and used mapping files')
                sys.exit()
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
    if decayRate < 0.01:
        traceability = 1
    else:
        traceability = 1 - ((decayPop * math.exp(decayRate * mlDist)) / (1 + decayPop * (math.exp(decayRate * mlDist) - 1)))
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


def generateNexusFile(nexus_file, taxonset, config):

    print_progress('Generating nexus file')

    with open(nexus_file, 'w') as fnew:
        fnew.write(f'#NEXUS\nbegin taxa;\n\tdimensions ntax={len(taxonset)}'
                   ';\n\ttaxlabels')
        for taxa in taxonset:
            fnew.write('\n\t' + taxa.replace(' ', '_'))
        fnew.write('\n;\nend;\n\nbegin trees;\n\ttree tree_1 =\n')
        f = open(config.species_tree).read()
        fnew.write(f)
        fnew.write('\nend;\n\nbegin figtree;\n')
        fnew.write('\tset appearance.backgroundColorAttribute="User Selection"'
                   ';\n\tset appearance.backgroundColour=#-1;\n\tset '
                   'appearance.branchColorAttribute="User Selection";\n\t'
                   'set appearance.branchLineWidth=3.0;\n\tset appearance.'
                   'foregroundColour=#-16777216;\n\tset appearance.'
                   'selectionColour=#-2144520576;\n\tset branchLabels.'
                   'colorAttribute="User Selection";\n\tset branchLabels.'
                   'displayAttribute="bootstrap";\n\tset branchLabels.'
                   'fontName="Times New Roman";\n\tset branchLabels.'
                   'fontSize=28;\n\tset branchLabels.fontStyle=1;\n\t'
                   'set branchLabels.isShown=true;\n\tset branchLabels.'
                   'significantDigits=4;\n\tset layout.expansion=0;\n\tset '
                   'layout.layoutType="POLAR";\n\tset layout.zoom=1100;\n\t'
                   'set nodeBars.barWidth=4.0;\n\tset nodeLabels.'
                   'colorAttribute="User Selection";\n\tset nodeLabels.'
                   'displayAttribute="Node ages";\n\tset nodeLabels.fontName='
                   '"sansserif";\n\tset nodeLabels.fontSize=14;\n\tset '
                   'nodeLabels.fontStyle=0;\n\tset nodeLabels.isShown=false;\n'
                   '\tset nodeLabels.significantDigits=4;\n\tset polarLayout.'
                   'alignTipLabels=false;\n\tset polarLayout.angularRange=0;'
                   '\n\tset polarLayout.rootAngle=0;\n\tset polarLayout.'
                   'rootLength=100;\n\tset polarLayout.showRoot=false;\n\tset '
                   'radialLayout.spread=0.0;\n\tset rectilinearLayout.'
                   'alignTipLabels=false;\n\tset rectilinearLayout.'
                   'curvature=0;\n\tset rectilinearLayout.rootLength=100;\n\t'
                   'set scale.offsetAge=0.0;\n\tset scale.rootAge=1.0;\n\t'
                   'set scale.scaleFactor=1.0;\n\tset scale.scaleRoot=false;'
                   '\n\tset scaleAxis.automaticScale=true;\n\tset scaleAxis.'
                   'fontSize=8.0;\n\tset scaleAxis.isShown=false;\n\tset '
                   'scaleAxis.lineWidth=1.0;\n\tset scaleAxis.majorTicks=1.0;'
                   '\n\tset scaleAxis.origin=0.0;\n\tset scaleAxis.'
                   'reverseAxis=false;\n\tset scaleAxis.showGrid=true;\n\tset '
                   'scaleAxis.significantDigits=4;\n\tset scaleBar.'
                   'automaticScale=true;\n\tset scaleBar.fontSize=10.0;\n\t'
                   'set scaleBar.isShown=true;\n\tset scaleBar.lineWidth=1.0;'
                   '\n\tset scaleBar.scaleRange=0.0;\n\tset scaleBar.'
                   'significantDigits=4;\n\tset tipLabels.colorAttribute='
                   '"User Selection";\n\tset tipLabels.displayAttribute='
                   '"Names";\n\tset tipLabels.fontName="Times New Roman";'
                   '\n\tset tipLabels.fontSize=18;\n\tset tipLabels.'
                   'fontStyle=1;\n\tset tipLabels.isShown=true;\n\tset '
                   'tipLabels.significantDigits=4;\n\tset trees.order=false;'
                   '\n\tset trees.orderType="increasing";\n\tset trees.'
                   'rooting=false;\n\tset trees.rootingType="User Selection";'
                   '\n\tset trees.transform=false;\n\tset trees.'
                   'transformType="cladogram";')
        fnew.write('\nend;')

