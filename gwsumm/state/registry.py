# -*- coding: utf-8 -*-
# Copyright (C) Duncan Macleod (2013)
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

"""Registry for GWSumm data states.

All Tabs should be registered for easy identification from the
configuration INI files
"""

from .. import (globalv, version)
from ..utils import re_quote

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version

__all__ = ['register_state', 'get_state']


def register_state(state, name=None, force=False):
    """Register a new `SummaryState` to the given ``name``

    Parameters
    ----------
    state : `SummaryState`
        defining Class for this state type.
    name : `str`, optional
        unique descriptive name for this type of state, must not
        contain any spaces, e.g. 'hveto'. If ``name=None``, the `Tab.type`
        class attribute of the given state will be used.
    force : `bool`
        overwrite existing registration for this type

    Raises
    ------
    ValueError
        if name is already registered and ``force`` not given as `True`
    """
    if name is None:
        name = state.name
    if not name in globalv.STATES or force:
        globalv.STATES[name] = state
    else:
        raise ValueError("Tab '%s' has already been registered to the %s "
                         "class" % (name, state.name))


def get_state(name):
    """Query the registry for the `SummaryState` registered to the given name
    """
    name = re_quote.sub('', name)
    try:
        return globalv.STATES[name]
    except KeyError:
        raise ValueError("No Tab registered with name '%s'" % name)
