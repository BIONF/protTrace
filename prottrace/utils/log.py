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

from time import (
    time,
    strftime,
    localtime
)


def print_progress(message, force_flush=False):
    print(f'###\t{message}\t###', flush=force_flush)


def print_warning(message):
    print(f'Warning: {message}')


def print_error(message):
    print(f'ERROR: {message}')


class time_report():
    __slots__ = ['start_time']

    def __init__(self):
        self.start_time = time()

    def time_passed(self):
        """ Returns the time passed since starting this instance in
        minutes. """
        return time() - self.start_time

    def format_secs(self):
        return f'{str(self.time_passed())} seconds'

    def mins_passed(self):
        """ Scales the time passed to minutes. """
        return self.time_passed() / 60

    def format_mins(self):
        """ Formats the time passed to show minutes. """
        return f'{str(self.mins_passed())} minutes'

    def hrs_passed(self):
        """ Scales the time passed to hours. """
        return self.time_passed() / 3600

    def format_hrs(self):
        """ Formats the time passed to show hours.. """
        return f'{str(self.hrs_passed())} hours'

    def format_passed(self, message, scale='sec'):
        """ Formats a given time with a specified message. """

        # The default time format is seconds.
        scaled_time = self.format_secs()
        if scale == 'min':
            scaled_time = self.format_mins()
        elif scale == 'hrs':
            scaled_time = self.format_hrs()

        return f'Time taken: {scaled_time} for {message}'

    def print_time(self, message, scale='sec'):
        """ Prints the time passed onto the user interface. """
        print_progress(self.format_passed(message, scale), force_flush=True)

    @staticmethod
    def format_verbose(message):
        """ Formats the current clock time verbosely. """
        formatted_time = strftime('%a, %d %b %Y %H:%M:%S +0000')
        print_progress(f'{formatted_time}: {localtime()}')
