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


def print_progress(message, force_flush=False):
    print('###\t{0}\t###'.format(message), flush=force_flush)


def print_warning(message):
    print('Warning: {0}'.format(message))


def print_error(message):
    print('ERROR: {0}'.format(message))


class time_report():
    __slots__ = ['start_time']

    def __init__(self):
        self.start_time = time.time()

    def time_passed(self):
        """ Returns the time passed since starting this instance in
        minutes. """
        return time.time() - self.start_time

    def format_secs(self):
        return '{0} seconds'.format(str(self.time_passed))

    def mins_passed(self):
        """ Scales the time passed to minutes. """
        return self.time_passed() / 60

    def format_mins(self):
        """ Formats the time passed to show minutes. """
        return '{0} minutes'.format(str(self.mins_passed()))

    def hrs_passed(self):
        """ Scales the time passed to hours. """
        return self.time_passed() / 3600

    def format_hrs(self):
        """ Formats the time passed to show hours.. """
        return '{0} hours'.format(str(self.hrs_passed()))

    def format_passed(self, message, scale='sec'):
        """ Formats a given time with a specified message. """

        # The default time format is seconds.
        scaled_time = self.format_secs()
        if scale == 'min':
            scaled_time = self.format_mins()
        elif scale == 'hrs':
            scaled_time = self.format_hrs()

        return '####\tTime taken: {0} for {1}\t####'.format(
            scaled_time, message)

    def print_time(self, message, scale='sec'):
        """ Prints the time passed onto the user interface. """
        print(self.format_passed(self, message, scale), flush=True)

    @staticmethod
    def format_verbose(message):
        """ Formats the current clock time verbosely. """
        print('#### {0}: {1} ####'.time.strftime(
            '%a, %d %b %Y %H:%M:%S +0000', time.localtime()))
