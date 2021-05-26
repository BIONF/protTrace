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

import sys
from pathlib import Path
import subprocess

# from species_distance.ml_matrix import main as ml_matrix_main

from utils.log import (
    print_progress,
    print_warning,
    print_error
)
from utils.file import generate_splitted_lines as gen_lines

# Module to reconstruct tree using RAxML and do scaling factor calculation


# Calculate the median of a set
def median(lst):
    even = (0 if len(lst) % 2 else 1) + 1
    half = int((len(lst) - 1) / 2)
    return sum(sorted(lst)[half:half + even]) / float(even)

# Rename the ortholog group sequences into oma ids
# (cut down every fasta header to 5 character long names)
# def rename_orth():
#   fnew = open('temp_orth_%s.fa' %protein_id, 'w')
#   for i in range(0, len(orth) - 1, 2):
#       fnew.write(orth[i][:6] + '\n' + orth[i+1] + '\n')
#   fnew.close()


# Perform MSA of the new renamed ortholog sequences file (MAFFT mafft)
def msa_convert(query_id, phy_file, config):
    if phy_file.exists():
        print_progress('Reusing existing alignment file: phy_file')
    else:
        print_progress('Generating MSA')
        subprocess.run([config.mafft, '--phylipout', f'ogSeqs_{query_id}.fa',
                        '>', f'ogSeqs_{query_id}.phy'], check=True, shell=True)
        # os.system('{0} --phylipout ogSeqs_{1}.fa > ogSeqs_{3}.phy'
        #          .format(mafft, protein_id, protein_id))


# added by ingo to get rid of raxml dependency
def run_iqtree(protein_id, reuse_cache, config):
    if (reuse_cache
            and Path(f'ogSeqs_{protein_id}.phy.treefile').exists()
            and Path(f'ogSeqs_{protein_id}.phy.ckp.gz').exists()):
        print_progress('ML tree already exists. Reusing it.')
    else:
        # os.system('rm -rf RAxML_*')
        cmd = [str(config.iqtree), '-nt', str(config.nr_processors), '-s',
               f'ogSeqs_{protein_id}.phy', '-m',
               config.aa_substitution_matrix, '-keep-ident', '-redo']
        if __debug__:
            print(cmd)
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError:
            print_error('IQTree has run into an issue for calculating the '
                        'substitution scaling factor')
            sys.exit()
        # os.system('{0} -nt {1} -s ogSeqs_{2}.phy -m {3} -keep-ident '
        #           '-redo'.format(config.iqtree, nr_proc, protein_id,
        # aaMatrix))


# Remove all the temp files generated
def rm_temp(protein_id):

    try:
        # The missing_ok parameter is added in version 3.8. Let's keep some
        # downwards compatibility.
        Path('temp_parameters_{protein_id}.txt').unlink()
        Path('maxLikDist_{protein_id}.txt').unlink()
    except FileNotFoundError:
        pass


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


def gen_orth_ml_dists(protein_id, orth):
    """ Reads the ML distance output file of IQTree and returns the
    pairwise species ML distances. """

    with Path(f'ogSeqs_{protein_id}.phy.mldist').open('r') as phy_dists:
        prot_dists = [line for line in gen_lines(phy_dists, sep=' ')]
        # The first line (list of tab-separated items, with 1 item) only
        # contains the number of species. We do make a sanity check first,
        # though.
        if (len(prot_dists[0]) == 1):
            if __debug__:
                print('DEBUG: The first line in the orthologous group species '
                      'distance file correctly contains the number of member '
                      'species.')
            del prot_dists[0]

    def get_pairwise_prot_dist(pairwise_distances, count_1, count_2):
        """ From a list of tab-separated lines, retrieve the float from
        position (count_1|count_2)."""

        # We do count from the 0th line in the file, but distances
        # start from the 0th + 1 column in each line.
        return float(pairwise_distances[count_2][count_1 + 1])

    # The order of columns is identical to the order of rows in respect to the
    # represented species.
    line_counts = range(len(prot_dists))
    for count_1 in line_counts:
        # The 0th column contains the species ID.
        spec1 = orth.member_prot_to_spec(prot_dists[count_1][0])

        # Generate the second iteration of the species pair list, without the
        # own species identifier. We do not care for same-species distances,
        # which is always 0. With the filter function, we spare the if-clause.
        for count_2 in filter(lambda count: count != count_1, line_counts):
            spec2 = orth.member_prot_to_spec(prot_dists[count_2][0])

            yield species_pair_distance(spec1, spec2,
                                        get_pairwise_prot_dist(prot_dists,
                                                               count_1,
                                                               count_2))


# Calculate the scaling factor based on maximum likelihood distances
def scaling_factor_max(protein_id, orth_group, cache_manager, config):

    # map_file = config.fdog_oma_tree_map
    scales = []
    # Generate maximum likelihood distance file for orthologs
    # protein_dist = Path(f'ogSeqs_{protein_id}.phy.mldist').read().split('\n')
    # ml_matrix_main(protein_dist, protein_id)

    print_progress('Collecting orthologous protein ML distances')

    # Collect the pairwise protein distance of all sequences in the
    # orthologous group. Each distance will be linked to a pair of species IDs.
    orth_dists = gen_orth_ml_dists(protein_id, orth_group)

    print_progress('Collecting pairwise species ML distances')

    # Collect the median pairwise species distance of all species represented
    # by protein sequences in this orthologous group.
    spec_dists = {pair.spec_set: pair for pair in
                  orth_group.gen_species_ml_dists(config)}

    print_progress('Drawing the median of all orthologous / species ML '
                   'distance ratios')

    try:
        scales = [orth.distance / spec_dists[orth.spec_set].distance
                  for orth in orth_dists]
    except KeyError:
        print_warning('There are not as many pairwise species distances '
                      'available as there are orthologous protein sequence '
                      'ML distances')
    # try:
    #     # orth_max = Path(f'maxLikDist_{protein_id}').read().split('\n')
    # species_max = config.species_maxLikMatrix.open('r').read().split('\n')
    #     # fdogFile = map_file.open('r').read().split('\n')
    #     for i in range(len(orth_max) - 1):
    #         line = orth_max[i].split('\t')[0]
    #         if '_' in line:
    #             species1 = line.split('_')[1]
    #         else:
    #             species1 = line

    #         #
    #         # This block is needed when max likelihood matrix do not have OMA
    #         # identifiers.
    #         #
    #         # print(species1)
    #         # for j in range(len(fdogFile) - 1):
    #         #	if species1 == fdogFile[j].split('\t')[3]:
    #         #		fdog1 = fdogFile[j].split('\t')[0]
    #         #		break

    #         for k in range(i + 1, len(orth_max) - 1):
    #             line = orth_max[k].split('\t')[0]
    #             if '_' in line:
    #                 species2 = line.split('_')[1]
    #             else:
    #                 species2 = line

    #             # ml_dist filenames
    #             ml_1 = cache_manager.cache_dir / f'{species1}_{species2}.lik'
    #             ml_2 = cache_manager.cache_dir / f'{species2}_{species1}.lik'

    #             # for j in range(len(fdogFile) - 1):
    #             #	if species2 == fdogFile[j].split('\t')[3]:
    #             #		fdog2 = fdogFile[j].split('\t')[0]
    #             #		break
    #             # print(species1, species2)

    #             maxDistOrth = float(orth_max[i].split('\t')[k + 1])
    #             maxDistSpecies = 0
    #             flag1 = True
    #             flag2 = True
    #             for species_max_i in range(len(species_max) - 1):
    #                 if species_max[species_max_i].split('\t')[0] == species1:
    #                     rowIndex = species_max_i
    #                     flag1 = False
    #           elif species_max[species_max_i].split('\t')[0] == species2:
    #                     columnIndex = species_max_i
    #                     flag2 = False
    #             # if not flag1 and not flag2:
    #             #	maxDistSpecies = float(species_max[rowIndex]
    #             #   .split('\t')[columnIndex])
    #             # else:
    #             #	maxDistSpecies = 1.0
    #             mlPresent = True
    #             if flag1 or flag2:
    #                 if cache_manager.in_cache(ml_1, True):
    #                     maxDistSpecies = float(ml_1.open('r')
    #                                            .read().split('\n')[0])
    #                 elif cache_manager.in_cache(ml_2, True):
    #                     maxDistSpecies = float(ml_2.open('r')
    #                                            .read().split('\n')[0])
    #                 else:
    #                     mlPresent = False
    #             else:
    #                 if (not species_max[rowIndex]
    #                         .split('\t')[columnIndex] == "NA"):
    #                     maxDistSpecies = float(species_max[rowIndex]
    #                                            .split('\t')[columnIndex])
    #                 else:
    #                     if cache_manager.in_cache(ml_1, True):
    #                         maxDistSpecies = float(ml_1.open('r')
    #                                                .read().split('\n')[0])
    #                     elif cache_manager.in_cache(ml_2, True):
    #                         maxDistSpecies = float(ml_2.open('r')
    #                                                .read().split('\n')[0])
    #                     else:
    #                         mlPresent = False

    #             if mlPresent and not maxDistSpecies == 0:
    #                 scales.append(maxDistOrth / maxDistSpecies)
    #             elif mlPresent and maxDistSpecies == 0:
    #                 pass
    # except Exception:
    #     print_error('Scaling factor calculation had an error')
    #     sys.exit('Maximum likelihood files are invalid!')

    if len(scales) >= 1:
        return median(scales)
    else:
        return None


# Main module for running tree reconstruction
def main(query, orthologs, evol_params, tree_file, cache, config):

    protein_id = query.id
    orth = orthologs
    phy_file = Path(f'ogSeqs_{protein_id}.phy')
    reuse_cache = config.reuse_cache
    del_tmp = config.delete_temp

    if orth.at_least_4_sequences():
        try:
            msa_convert(protein_id, phy_file, config)
            run_iqtree(protein_id, reuse_cache, config)
            if reuse_cache and evol_params.file_exists():
                print_progress('Scaling factor already exists. Reusing it.')
            else:
                print_progress('Calculating the substitution scaling factor.')
                evol_params.set_scaling_factor(scaling_factor_max(protein_id,
                                                                  orth,
                                                                  cache,
                                                                  config))
            if del_tmp:
                rm_temp(protein_id)
        except Exception as e:
            raise e
            print_error('Some step in the tree reconstruction was invalid!')
            pass
    else:
        print_progress('Using default scaling factor: '
                       f'{str(evol_params.scaling_factor)}')
