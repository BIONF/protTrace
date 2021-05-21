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

from pathlib import Path
from utils.log import print_warning

# This module defines classes for handling evolutionary model parameters.


class evol_parameters:
    """ A container of evolutionary model parameters for ProtTrace. Parameters
    are first initialized using thir default values from the ProtTrace config.
    Later insertions then turn the default flag off. """
    __slots__ = ['_query_id', '_indel_rate', '_indel_length_distribution',
                 '_scaling_factor']

    def __init__(self, query, config):
        self._query_id = query.id

        self._indel_rate = evol_parameter(
            config.default_indel)
        self._indel_length_distribution = evol_parameter(
            config.default_indel_distribution)
        self._scaling_factor = evol_parameter(
            config.default_scaling_factor)

    @property
    def indel_rate(self):
        return self._indel_rate.value

    @property
    def indel_length_distribution(self):
        return self._indel_length_distribution.value

    @property
    def scaling_factor(self):
        return self._scaling_factor.value

    def set_indel_rate(self, rate):
        if not check_none('indel rate', rate):
            self._indel_rate.set_value(rate)

    def set_indel_length_distribution(self, dist):
        if not check_none('indel length distribution', dist):
            self._indel_length_distribution.set_value(dist)

    def set_scaling_factor(self, scale):
        if not check_none('substitution scaling factor', scale):
            self._scaling_factor.set_value(scale)

    def file_exists(self):
        return Path(self._query_id + '_parameters.txt').exists()

    def write_to_file(self):
        filepath = Path(self._query_id + '_parameters.txt')
        content = (''
                   f'indel rate:\t{str(self.indel_rate)}\n'
                   'indel length distribution:\t'
                   f'{str(self.indel_length_distribution)}\n'
                   'substitution scaling factor:\t'
                   f'{str(self.scaling_factor)}\n'
                   )
        with filepath.open('w') as evol_file:
            evol_file.write(content)

        return filepath


def check_none(name, value):
    if value is None:
        print_warning(f'The {name} was invalid. Reverting to default')
    return value is None


class evol_parameter:
    """ Contains the value of an evolutionary parameter. Is expected to be
    initialized with default values from the configuration. Later insertions
    can use set_value to also turn the default flag to False. """
    __slots__ = ['value', 'is_default']

    def __init__(self, value):
        self.value = value
        self.is_default = True

    def set_value(self, value):
        self.value = value
        self.is_default = False

    def __str__(self):
        return str(self.value)
