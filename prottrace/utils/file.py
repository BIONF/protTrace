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

""" This python module contains small code snippets to process files. """

import os
from pathlib import Path


class dir_move():
    __slots__ = ['current', 'previous']

    def __init__(self, target):
        self.previous = Path.cwd()
        self.current = target
        if not self.current.exists():
            self.current.mkdir()
        os.chdir(target)

    def reverse(self):
        os.chdir(self.previous)


def generate_splitted_lines_with_pos(handler):
    while True:
        previous_pos = handler.tell()
        line = handler.readline()
        if not line:
            break
        yield (previous_pos, line.rstrip().split('\t'))


def generate_splitted_lines(lines, comment='#'):
    for line in lines:
        if comment not in line:
            yield line.rstrip().split('\t')
