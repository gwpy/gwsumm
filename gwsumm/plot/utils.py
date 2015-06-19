# -*- coding: utf-8 -*-
# Copyright (C) Duncan Macleod (2013-2015)
#
# This file is part of GWSumm.
#
# GWSumm is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GWSumm is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GWSumm.  If not, see <http://www.gnu.org/licenses/>.

"""Utilies for GWSumm plotting
"""

from gwpy.plotter.table import get_column_string

from .. import version

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version

COLUMN_LABEL = {
    'peal_frequency': r"Frequency [Hz]",
    'central_freq': r"Frequency [Hz]",
    'frequency': r"Frequency [Hz]",
    'mchirp': r"Chirp mass [M$_\odot$]",
    'new_snr': r"$\chi^2$-weighted signal-to-noise ratio (New SNR)",
    'peak_frequency': r"Frequency [Hz]",
    'rho': r"$\rho$",
    'snr': r"Signal-to-noise ratio (SNR)",
    'template_duration': r"Template duration [s]",
}


def get_column_label(column):
    try:
        return COLUMN_LABEL.get(column)
    except KeyError:
        return get_column_string(column)
