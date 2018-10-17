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

"""Registry for `states <SummaryState>`.
"""

from .. import globalv
from ..utils import re_quote

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'

__all__ = ['register_state', 'get_state', 'get_states']


def register_state(state, key=None, force=False):
    """Register a new `SummaryState` to the given ``key``

    Parameters
    ----------
    state : `SummaryState`
        defining Class for this state type.
    key : `str`, optional
        unique descriptive name for the `SummaryState` to be registered.
        If ``key=None``, the :attr:`~SummaryState.key`
        attribute of the given state will be used.
    force : `bool`
        overwrite existing registration for this key

    Raises
    ------
    ValueError
        if key is already registered and ``force`` not given as `True`
    """
    if key is None:
        key = state.key
    key = key.lower()
    if key not in globalv.STATES or force:
        globalv.STATES[key] = state
        return state
    raise ValueError("State %r has already been registered." % key)


def get_state(key):
    """Query the registry for the `SummaryState` registered to the given key

    Parameters
    ----------
    key : `str`
        registered key of desired `SummaryState`. This may not match the
        `~SummaryState.name` attribute` if the state was registered with
        a different key.

    Returns
    -------
    state : `SummaryState`
        the `SummaryState` registered with the given key

    Raises
    ------
    ValueError:
        if the ``key`` doesn't map to a registered `SummaryState`
    """
    key = re_quote.sub('', key)
    try:
        return globalv.STATES[key.lower()]
    except KeyError:
        raise ValueError("No SummaryState registered with name '%s'" % key)


def get_states(keys=set()):
    """Query the registry for a list of states (defaults to all)

    Parameters
    ----------
    keys : `set` of `str`
        the set of state keys to query in the registry

    Returns
    -------
    states : `dict`
        a `dict` of (``key``, `SummaryState`) pairs

    Raises
    ------
    ValueError:
        if any of the ``keys`` doesn't map to a registered `SummaryState`
    """
    if not keys:
        return globalv.STATES.copy()
    else:
        return dict((key, get_state(key)) for key in keys)
