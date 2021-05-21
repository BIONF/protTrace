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
from utils.log import print_warning


# Script to transform MSA into indel blocks
# The output file is used as input for IQTree24 to calculate indel rates
def preprocess_phylip_file(phy_file):
    with phy_file.open('r') as phy_f:
        phy_file_content = phy_f.read().split('\n')
    totalSpecies = int(phy_file_content[0].split()[0])
    phy_edit = []
    c = 0
    for line in phy_file_content[1:totalSpecies + 1]:
        complete_sequence = line[11:].replace(' ', '')
        addLine = c + totalSpecies + 2
        while addLine < len(phy_file_content):
            complete_sequence += phy_file_content[addLine][11:].replace(' ',
                                                                        '')
            addLine += totalSpecies + 1
        phy_edit.append(complete_sequence)
        c += 1

    return phy_edit


# Calculates the number of indel events occuring in the alignment
# Divides the complete alignment based on indel blocks and stores the column
# number where indel occurred
def calculate_indel_blocks(f):
    indelBlocks = []

    indelRows = []
    for i in range(len(f[0])):
        indelStart = False
        for j in range(len(f)):
            if f[j][i] == '-':
                if j not in indelRows:
                    indelRows.append(j)
                    indelStart = True
            else:
                if j in indelRows:
                    indelRows.remove(j)
        if indelStart:
            indelBlocks.append(i)

    return indelBlocks


# Creates the transformed MSA file
def create_transform_align(trans_file, phy_file, edit_phy_file, indels_pos):
    with trans_file.open('w') as fnew:
        with phy_file.open('r') as phy_f:
            phy = phy_f.read().split('\n')
        totalSpecies = int(phy[0].split()[0])
        fnew.write(str(totalSpecies) + ' ' + str(len(indels_pos)) + '\n')
        c = 0

        def add_gap_count(handle, edit_file, c, pos_from, pos_to):
            handle.write(' ' +
                         str(edit_file[c][pos_from:pos_to].count('-')))

        def add_total_gap_count(handle, edit_file, c):
            handle.write(' ' + str(edit_file[c].count('-')))

        for line in phy[1:totalSpecies + 1]:
            fnew.write(line.split()[0])
            if len(indels_pos) == 1:
                start_pos = int(indels_pos[0])
                add_gap_count(fnew, edit_phy_file, c, start_pos,
                              len(edit_phy_file[c]))
            elif len(indels_pos) > 1:
                for i in range(len(indels_pos) - 1):
                    start_pos = int(indels_pos[i])
                    stop_pos = int(indels_pos[i + 1])
                    add_gap_count(fnew, edit_phy_file, c, start_pos, stop_pos)
                add_gap_count(fnew, edit_phy_file, c, stop_pos,
                              len(edit_phy_file[c]))
            else:
                add_total_gap_count(fnew, edit_phy_file, c)
            fnew.write('\n')
            c += 1


def post_transform_align(trans):
    t = open(trans).read()
    t0 = open(trans).read().split('\n')

    result = open(trans, 'w')
    result.write(t0[0] + '\n')
    try:
        stateOld = []

        for i in range(1, len(t0) - 1):
            statesList = t0[i].split()[1:]
            for s in statesList:
                if not int(s) in stateOld:
                    stateOld.append(int(s))

        stateOld.sort()
        stateOld = stateOld[::-1]

        for i in range(len(stateOld)):
            # DEBUG
            # print(len(stateOld) - i, stateOld[i])

            # This is the core operation to convert existing data into relevant
            # indel information.
            t = t.replace(' ' + str(stateOld[i]) + ' ',
                          ' ' + str(len(stateOld) - 1 - i) + ' ').replace(
                              ' ' + str(stateOld[i]) + '\n',
                              ' ' + str(len(stateOld) - 1 - i) + '\n')
    except Exception:
        print_warning('Error while editing the transformed alignment file')
        pass

    for i in range(1, len(t.split('\n')) - 1):
        result.write(t.split('\n')[i] + '\n')


def main(phy_file, trans_file):
    if not os.stat(phy_file).st_size == 0:
        edit_file = preprocess_phylip_file(phy_file)
        indel_blocks_pos = calculate_indel_blocks(edit_file)
        create_transform_align(trans_file, phy_file, edit_file,
                               indel_blocks_pos)
        post_transform_align(trans_file)
        return len(indel_blocks_pos)
    else:
        print_warning('Transformed alignment cannot be created for empty '
                      'phylip files')
