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

"""Registry for GWSumm data tabs.

All Tabs should be registered for easy identification from the
configuration INI files
"""

from ..utils import re_quote

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'

__all__ = ['register_tab', 'get_tab']

_TABS = {}


def register_tab(tab, name=None, force=False):
    """Register a new summary `Tab` to the given ``name``

    Parameters
    ----------
    tab : `type`
        defining Class for this tab type.
    name : `str`, optional
        unique descriptive name for this type of tab, must not
        contain any spaces, e.g. 'hveto'. If ``name=None``, the `Tab.type`
        class attribute of the given tab will be used.
    force : `bool`
        overwrite existing registration for this type

    Raises
    ------
    ValueError
        if name is already registered and ``force`` not given as `True`
    """
    if name is None:
        name = tab.type
    if name not in _TABS or force:
        _TABS[name] = tab
    else:
        raise ValueError("Tab '%s' has already been registered to the %s "
                         "class" % (name, _TABS[name].__name__))


def get_tab(name):
    """Query the registry for the tab class registered to the given
    name
    """
    name = re_quote.sub('', name)
    try:
        return _TABS[name]
    except KeyError:
        raise ValueError("No Tab registered with name '%s'" % name)
