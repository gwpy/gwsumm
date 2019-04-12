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

"""Thin wrapper around configparser
"""

import configparser
import os.path
import re
from collections import OrderedDict
from importlib import import_module

from six import string_types
from six.moves import (
    StringIO,
    http_client as httplib,
)

# import these for evaluating lambda expressions in the configuration file
import math  # noqa: F401
import numpy  # noqa: F401

from matplotlib import rcParams

from astropy import units

from gwpy.detector import (Channel, ChannelList)
from gwpy.time import tconvert

from .html import (get_css, get_js)
from .utils import (nat_sorted, re_cchar, re_quote, safe_eval, OBSERVATORY_MAP)
from .channels import (get_channels, split as split_channels,
                       update_channel_params)

__all__ = [
    'GWSummConfigParser',
]


class GWSummConfigParser(configparser.ConfigParser):
    # preserve case in options
    optionxform = str
    # disable colon separator
    OPTCRE = re.compile(
        r'(?P<option>[^=\s][^=]*)\s*(?P<vi>[=])\s*(?P<value>.*)$')
    # set interpolation match
    _interpvar_re = re.compile(r"%\(([^)]+)\)s")

    def __init__(self, *args, **kwargs):
        kwargs.setdefault('dict_type', OrderedDict)
        configparser.ConfigParser.__init__(self, *args, **kwargs)
    __init__.__doc__ = configparser.ConfigParser.__init__.__doc__

    def read_file(self, *args, **kwargs):
        try:
            return configparser.ConfigParser.read_file(self, *args, **kwargs)
        except AttributeError:  # python < 3
            return self.readfp(*args, **kwargs)

    def ndoptions(self, section, **kwargs):
        options = configparser.ConfigParser.options(self, section, **kwargs)
        return [o for o in options if o not in self._defaults]

    def nditems(self, section, **kwargs):
        items = configparser.ConfigParser.items(self, section, **kwargs)
        return [i for i in items if i[0] not in self._defaults]

    def read(self, filenames):
        readok = configparser.ConfigParser.read(self, filenames)
        if isinstance(filenames, string_types):
            filenames = filenames.split(',')
        for fp in filenames:
            if fp not in readok:
                raise IOError("Cannot read file: %s" % fp)
        self.files = list(map(os.path.abspath, readok))
        return readok

    @classmethod
    def from_configparser(cls, cp):
        """Copy an existing :class:`~configparser.ConfigParser`.
        """
        # if its already what we need, just return this instance
        if isinstance(cp, cls):
            return cp
        # set up temporary buffer
        buf = StringIO()
        # write to buffer
        cp.write(buf)
        buf.seek(0)
        # read new GWSummConfigParser
        new = cls()
        new.read_file(buf)
        return new

    def __repr__(self):
        return '<GWSummConfigParser()>'

    def interpolate_section_names(self, **kwargs):
        """Interpolate a specific key in a section name using val
        """
        for section in self.sections():
            s = section
            for key in self._interpvar_re.findall(section):
                try:
                    val = kwargs[key]
                except KeyError:
                    raise configparser.InterpolationMissingOptionError(
                        '[%s]' % section, section, key, '%%(%s)s' % key)
                s = section % {key: val}
            self._sections[s] = self._sections.pop(section)

    def set_ifo_options(self, ifo, observatory=None,
                        section=configparser.DEFAULTSECT):
        """Set configurations options in [DEFAULT] based on the given `ifo`

        The following options are set

        - `IFO` - the two-character interferometer prefix, e.g. ``L1``
        - `ifo` - the two-character interferometer prefix in lower-case,
           e.g. ``l1``
        - `SITE` - the single-character site ID, e.g. ``L``
        - `site` - the single-character site ID,n lower-case e.g. ``l``

        Additionally, if `observatory` is given, or the `ifo` matches known
        observatories, the following option is set

        - `observatory` - the name of the observatory, e.g. ``LIGO Livingston``

        """
        if observatory is None:
            observatory = OBSERVATORY_MAP.get(ifo)
        self.set(section, 'IFO', ifo)
        self.set(section, 'ifo', ifo.lower())
        self.set(section, 'SITE', ifo[0].upper())
        self.set(section, 'site', ifo[0].lower())
        if observatory is not None:
            self.set(section, 'observatory', observatory)

    def set_date_options(self, start, end, section=configparser.DEFAULTSECT):
        """Set datetime options in [DEFAULT] based on the given times

        The following options are set

        - `gps-start-time` - the integer GPS start time of this job
        - `gps-end-time` - the integer GPS end time of this job
        - `yyyy` - the four-digit year of the start date
        - `mm` - the two-digit month of the start date
        - `dd` - the two-digit day-of-month of the start date
        - `yyyymm` - the six-digit year and month of the start date
        - `yyyymmdd` - the eight-digit year-month-day of the start date
        - `duration` - the duration of the job (seconds)

        Additionally, if LAL is available, the following extra options
        are also set

        - `leap-seconds` - the number of leap seconds for the start date
        - `gps-start-time-noleap` - the leap-corrected integer GPS start time
        - `gps-end-time-noleap` - the leap-corrected integer GPS end time

        """
        utc = tconvert(start)
        self.set(section, 'gps-start-time', str(int(start)))
        self.set(section, 'gps-end-time', str(int(end)))
        self.set(section, 'yyyy', utc.strftime('%Y'))
        self.set(section, 'yy', utc.strftime('%y'))
        self.set(section, 'mm', utc.strftime('%m'))
        self.set(section, 'dd', utc.strftime('%d'))
        self.set(section, 'yyyymm', utc.strftime('%Y%m'))
        self.set(section, 'yyyymmdd', utc.strftime('%Y%m%d'))
        self.set(section, 'duration', str(int(end - start)))
        try:
            from lal import GPSLeapSeconds as leap_seconds
        except ImportError:
            pass
        else:
            nleap = leap_seconds(int(start))
            self.set(section, 'leap-seconds', str(nleap))
            self.set(section, 'gps-start-time-noleap', str(int(start) - nleap))
            self.set(section, 'gps-end-time-noleap', str(int(end) - nleap))

    def finalize(self):
        """Finalize this `GWSummConfigParser` by running all the loaders

        This method is just a shortcut to run each of the following

        .. autosummary::

           ~GWSummConfigParser.load_plugins
           ~GWSummConfigParser.load_units
           ~GWSummConfigParser.load_channels
           ~GWSummConfigParser.load_states

        """
        self.load_plugins()
        self.load_units()
        self.load_channels()
        self.load_states()

    def load_plugins(self):
        """Load all plugin modules as defined in the [plugins] section
        """
        mods = []
        try:
            plugins = self.ndoptions('plugins')
        except configparser.NoSectionError:
            pass
        else:
            for plugin in plugins:
                mods.append(import_module(plugin))
        return mods

    def load_units(self):
        """Load all unit definitions as defined in the [units] section
        """
        try:
            customunits = self.nditems('units')
        except configparser.NoSectionError:
            return []
        else:
            new_ = []
            for unit, b in customunits:
                if b.lower() == 'dimensionless':
                    b = ''
                new_.append(units.def_unit([unit], units.Unit(b)))
            units.add_enabled_units(new_)
            return new_

    def load_channels(self):
        """Load all channel definitions as given in the selfuration

        Channels are loaded from sections named [channels-...] or
        those sections whose name is a channel name in itself
        """
        channels_re = re.compile(r'channels[-\s]')
        # parse channel grouns into individual sections
        for section in filter(channels_re.match, self.sections()):
            names = split_channels(self.get(section, 'channels'))
            items = dict(self.nditems(section, raw=True))
            items.pop('channels')
            for name in names:
                name = name.strip(' \n')
                if not self.has_section(name):
                    self.add_section(name)
                for key, val in items.items():
                    if not self.has_option(name, key):
                        self.set(name, key, val)

        # read all channels
        raw = set()
        trend = set()
        for section in self.sections():
            try:
                m = Channel.MATCH.match(section).groupdict()
            except AttributeError:
                continue
            else:
                if not m['ifo']:
                    continue
            if m['trend']:
                trend.add(section)
            else:
                raw.add(section)

        channels = ChannelList()
        for group in [raw, trend]:
            try:
                newchannels = get_channels(group)
            except httplib.HTTPException:
                newchannels = []

            # read custom channel definitions
            for channel, section in zip(newchannels, group):
                for key, val in nat_sorted(self.nditems(section),
                                           key=lambda x: x[0]):
                    key = re_cchar.sub('_', key.rstrip())
                    if key.isdigit():
                        if not hasattr(channel, 'bits'):
                            channel.bits = []
                        while len(channel.bits) < int(key):
                            channel.bits.append(None)
                        if val.startswith('r"') or val.startswith('r\''):
                            val = eval(val)
                        channel.bits.append(str(val))
                    else:
                        setattr(channel, key, safe_eval(val.rstrip()))
                channels.append(channel)

        # rebuild channel list with new parameters
        update_channel_params()

        return channels

    def load_states(self, section='states'):
        """Read and format a list of `SummaryState` definitions from the
        given :class:`~configparser.ConfigParser`
        """
        from .state import (register_state, SummaryState,
                            ALLSTATE, generate_all_state, get_state)
        # parse the [states] section into individual state definitions
        try:
            states = dict(self.nditems(section))
        except configparser.NoSectionError:
            self.add_section(section)
            states = {}
        for state in states:
            if not (self.has_section('state-%s' % state) or
                    self.has_section('state %s' % state)):
                section = 'state-%s' % state
                self.add_section(section)
                self.set(section, 'name', state)
                self.set(section, 'definition', states[state])

        # parse each state section into a new state
        states = []
        for section in self.sections():
            if re.match(r'state[-\s]', section):
                states.append(register_state(
                    SummaryState.from_ini(self, section)))

        # register All state
        start = self.getint(section, 'gps-start-time')
        end = self.getint(section, 'gps-end-time')
        try:
            all_ = generate_all_state(start, end)
        except ValueError:  # user defined an all state themselves
            all_ = get_state(ALLSTATE)
        else:
            states.insert(0, all_)

        return states

    def load_rcParams(self, section='rcParams'):
        """Load custom `matplotlib.rcParams` for plots in this analysis
        """
        try:
            new = dict(self.nditems(section))
        except configparser.NoSectionError:
            return dict()
        else:
            for key, value in new.items():
                new[key] = safe_eval(value)
            rcParams.update(new)
            return new

    def get_css(self, section='html'):
        # get critical CSS
        css = get_css().copy()
        for key in css:
            try:
                css[key] = self.get(section, '%s-css' % key)
            except configparser.NoSectionError:  # no overrides are present
                return list(css.values())
            except configparser.NoOptionError:
                continue
        files = list(css.values())
        # get extra CSS
        try:
            extras = self.get(section, 'extra-css')
        except configparser.NoOptionError:
            return files
        else:
            files.extend(map(lambda x: re_quote.sub('', x), extras.split(',')))
        return files

    def get_javascript(self, section='html'):
        # get critical JS
        js = get_js().copy()
        for key in js:
            try:
                js[key] = self.get(section, '%s-js' % key)
            except configparser.NoSectionError:
                return list(js.values())
            except configparser.NoOptionError:
                continue
        files = list(js.values())
        # get extra CSS
        try:
            extras = self.get(section, 'extra-js')
        except configparser.NoOptionError:
            return files
        else:
            files.extend(map(lambda x: re_quote.sub('', x), extras.split(',')))
        return files
