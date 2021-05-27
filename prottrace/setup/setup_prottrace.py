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
import os

from argparse import ArgumentParser
from subprocess import run
from pathlib import Path


def create_directories(*dirs):
    """ Tests whether each path exists and creates the directory if it
        does not. """
    for path in dirs:
        if not path.is_dir():
            os.mkdir(path)


def main():

    def load_arguments():
        version = '1.8.4'
        parser = ArgumentParser(description='The setup for ProtTrace '
                                         'version {0}. Information can be '
                                         'changed later in the configuration '
                                         'file.'.format(version))

        meta_arguments = parser.add_argument_group('General arguments')
        meta_arguments.add_argument('--conda', action='store_true',
                                    default=False,
                                    help='Setup ProtTrace within a conda '
                                    'environment.')
        meta_arguments.add_argument('-p', '--processors', type=int, nargs='?',
                                    default=os.cpu_count(),
                                    help='Restrict the number of available '
                                    'processors.')

        orthologs = parser.add_argument_group('Orthologous sequences input')
        orthologs.add_argument('--fdog', action='store_true',
                                default=False,
                                help='Use HaMStR for detecting orthologs.')
        output = parser.add_argument_group('Output arguments')
        output.add_argument('-c', '--configout', type=Path,
                             default=Path('./prog.config'),
                             help='The configuration file path.')
        output.add_argument('-o', '--out', type=Path,
                             default=Path('./output/'),
                             help='The output directory.')
        output.add_argument('--cache', type=Path, default=Path('./cache/'),
                             help='The cache directory.')
        output.add_argument('--dist_dir', type=Path,
                             default=Path('./distances/'),
                             help='The temporary pairwise distance '
                             'computation path.')

        return parser.parse_args()

    args = load_arguments()
    conda = args.conda
    nr_procs = args.processors
    config_path = args.configout

    create_directories(args.out, args.cache, args.dist_dir)

    setup_path = Path('./bin/setup/setup.sh')
    run([str(setup_path)])

#    def organize_config(config_path, args):
#        """ Fills the ProtTrace configuration file with the determined
#            arguments. """
#

if __name__ == '__main__':
    main()
