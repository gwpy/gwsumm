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

import hashlib
import itertools
import re

from matplotlib import rcParams

from gwpy.plot.tex import label_to_latex
from gwpy.plot.utils import (  # noqa: F401
    FIGURE_PARAMS,
    AXES_PARAMS,
)

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'

# -- plotting parameters ------------------------------------------------------

LINE_PARAMS = [
    'linewidth', 'linestyle', 'color', 'label', 'alpha', 'rasterized',
    'zorder',
]
COLLECTION_PARAMS = [
    'cmap', 'vmin', 'vmax', 'marker', 's', 'norm', 'rasterized',
]
IMAGE_PARAMS = [
    'imshow', 'cmap', 'vmin', 'vmax', 'norm', 'rasterized', 'extent',
    'origin', 'interpolation', 'aspect',
]
HIST_PARAMS = [
    'bins', 'range', 'normed', 'weights', 'cumulative', 'bottom',
    'histtype', 'align', 'orientation', 'rwidth', 'log', 'color',
    'label', 'stacked', 'logbins',
]
LEGEND_PARAMS = [
    'loc', 'borderaxespad', 'ncol',
]
ARTIST_PARAMS = set(itertools.chain.from_iterable([
    LINE_PARAMS,
    COLLECTION_PARAMS,
    IMAGE_PARAMS,
    HIST_PARAMS,
]))

# -- default labels for table columns -----------------------------------------

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
        return COLUMN_LABEL[column]
    except KeyError:
        return get_column_string(column)


def get_column_string(column):
    r"""
    Format the string columnName (e.g. xml table column) into latex format for
    an axis label.

    Parameters
    ----------
    column : `str`
        string to format

    Examples
    --------
    >>> get_column_string('snr')
    'SNR'
    >>> get_column_string('bank_chisq_dof')
    r'Bank $\chi^2$ DOF'
    """
    acro = ['snr', 'ra', 'dof', 'id', 'ms', 'far']
    greek = ['alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'eta',
             'theta', 'iota', 'kappa', 'lamda', 'mu', 'nu', 'xi', 'omicron',
             'pi', 'rho', 'sigma', 'tau', 'upsilon', 'phi', 'chi', 'psi',
             'omega']
    unit = ['ns']
    sub = ['flow', 'fhigh', 'hrss', 'mtotal', 'mchirp']

    tex = rcParams['text.usetex']

    words = []
    for word in re.split(r'\s', column):
        if word.isupper():
            words.append(word)
        else:
            words.extend(re.split(r'_', word))

    # parse words
    for i, word in enumerate(words):
        # get acronym in lower case
        if word in acro:
            words[i] = word.upper()
        # get numerical unit
        elif word in unit:
            words[i] = '$(%s)$' % word
        # get character with subscript text
        elif word in sub and tex:
            words[i] = r'%s$_{\mbox{\small %s}}$' % (word[0], word[1:])
        # get greek word
        elif word in greek and tex:
            words[i] = r'$\%s$' % word
        # get starting with greek word
        elif re.match(r'(%s)' % '|'.join(greek), word) and tex:
            if word[-1].isdigit():
                words[i] = r'$\%s_{%s}$''' % tuple(
                    re.findall(r"[a-zA-Z]+|\d+", word))
            elif word.endswith('sq'):
                words[i] = r'$\%s^2$' % word.rstrip('sq')
        # get everything else
        else:
            if word.isupper():
                words[i] = word
            else:
                words[i] = word.title()
            # escape underscore
            words[i] = usetex_tex(re.sub(r'(?<!\\)_', r'\_', words[i]))
    return ' '.join(words)


def usetex_tex(text):
    """Format text for TeX if `text.usetex` is True
    """
    if rcParams['text.usetex']:
        return label_to_latex(text)
    return text


def hash(string, num=6):
    """Generate an N-character hash string based using string to initialise

    Parameters
    ----------
    string : `str`
        the initialisation string

    num : `int`, optional
        the length of the hash to produce

    Returns
    -------
    hash : `str`
        the new hash

    Examples
    --------
    >>> from gwsumm.plot.utils import hash
    >>> print(hash("I love gravitational waves"))
    80c897
    """
    return hashlib.md5(string.encode("utf-8")).hexdigest()[:num]
