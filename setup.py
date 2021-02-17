#!/usr/bin/env python
# -*- coding: utf-8 -*-
#######################################################################
#  Copyright (C) 2021 Arpit Jain, Dominik Perisa,
#  Prof. Dr. Ingo Ebersberger <ebersberger@bio.uni-frankfurt.de>
#
#  ProtTrace is a free software: you can redistribute it and/or
#  modify it under the terms of the GNU General Public License as
#  published by the Free Software Foundation, either version 3 of
#  the License, or (at your option) any later version.
#
#  ProtTrace is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with ProtTrace. If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    long_description=long_description,
    long_description_content_type="text/markdown",
)
