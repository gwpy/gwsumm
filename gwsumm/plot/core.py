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


"""Parse, define, and geenrate plots as requests through the configuration
for GWSumm
"""

import hashlib
import os.path
import re
import warnings
from math import (floor, ceil)
from urlparse import urlparse

try:
    from collections import OrderedDict
except ImportError:
    from astropy.utils import OrderedDict

from gwpy.segments import Segment
from gwpy.detector import ChannelList
from gwpy.plotter.utils import rUNDERSCORE

from .registry import register_plot
from .. import globalv
from ..data import get_channel
from ..utils import split_channels

__all__ = ['SummaryPlot', 'DataPlot']

re_cchar = re.compile("[\W\s_]+")


class SummaryPlot(object):
    """An image to displayed in GWSumm HTML output.

    Parameters
    ----------
    href : `str`, optional
        The IMG URL for this `SummaryPlot`.
    new : `bool`, optional
        `bool` flag whether this is a new plot to be processed (`True`),
        of that the output already exists on disk (`False`).

    Notes
    -----
    This `class` is a stub, designed to make creating detailed `SummaryPlot`
    classes easier.
    """
    type = None
    _threadsafe = True

    def __init__(self, href=None, new=True):
        self.href = href
        self.new = new

    @property
    def href(self):
        """HTML <img> href attribute for this `SummaryPlot`.
        """
        return self._href

    @href.setter
    def href(self, url):
        if url is None:
            self._href = None
        elif urlparse(url).netloc:
            self._href = url
        else:
            self._href = os.path.normpath(url)

    @property
    def new(self):
        """Flag whether this is a new plot or, already exists.

        Set new=False to skip actually processing this `SummaryPlot`, and
        just link to the outputfile.
        """
        return self._new

    @new.setter
    def new(self, isnew):
        self._new = bool(isnew)

    # ------------------------------------------------------------------------
    # TabSummaryPlot methods

    @classmethod
    def from_ini(cls, *args, **kwargs):
        """Define a new `SummaryPlot` from a an INI-format `ConfigParser`
        section.
        """
        raise NotImplementedError("Sub-classes should provide this method.")

    def __eq__(self, other):
        """Compare this `SummaryPlot` to another.

        Returns
        -------
        `True`
            if the metadata for ``self`` match those of ``other``
        `False`
            otherwise
        """
        if not isinstance(other, self.__class__):
            return False
        if not self.href == other.href:
            return False
        return True

    def __repr__(self):
        return '<%s(%s)>' % (self.__class__.__name__, self.href)

    def __str__(self):
        return str(self.href)

register_plot(SummaryPlot)


class DataPlot(SummaryPlot):
    """A `SummaryPlot` from instrumental data.

    Parameters
    ----------
    channels : `list`
        a list of channel names that define the data sources for this
        `DataPlot`
    start : `float`
        GPS start time of this `DataPlot`.
    end : `float`
        GPS end time of this `DataPlot`.
    tag : `str`
        a descriptive tag for this `TabSummaryPlot`, used as part of the output
        file name
    outdir : `str`
        output directory path for this `TabSummaryPlot`, defaults to the current
        directory
    href : `str`
        custom URL for this plot to link towards.
    **kwargs
        all other keyword arguments to be passed to this plot's
        :meth:`process` method.

    Notes
    -----
    All sub-classes of this object must provide the following methods

    =======================  ==============================================
    :meth:`add_data_source`  routine for appending data sources to the plot
    =======================  ==============================================
    """
    #: name for TabSummaryPlot subclass
    type = 'data'
    #: dict of default plotting kwargs
    defaults = {}

    def __init__(self, channels, start, end, state=None, outdir='.',
                 tag=None, href=None, new=True, all_data=False, **pargs):
        super(DataPlot, self).__init__(href=href, new=new)
        self.channels = channels
        self.span = (start, end)
        self.state = state
        self._outdir = outdir
        if tag is not None:
            self.tag = tag
        self.all_data = all_data
        self.pargs = self.defaults.copy()
        self.pargs.update(pargs)
        self.plot = None

    # ------------------------------------------------------------------------
    # TabSummaryPlot properties

    @property
    def span(self):
        """The GPS [start, stop) interval for this `TabSummaryPlot`.
        """
        return self._span

    @span.setter
    def span(self, seg):
        self._span = Segment(*seg)

    @property
    def start(self):
        return self.span[0]

    @property
    def end(self):
        return self.span[1]

    @property
    def state(self):
        """`~gwsumm.state.SummaryState` defining validity of this
        `DataPlot`.
        """
        return self._state

    @state.setter
    def state(self, state_):
        if isinstance(state_, (unicode, str)):
            self._state = globalv.STATES[state_]
        else:
            self._state = state_

    @property
    def channels(self):
        """List of data-source
        :class:`Channels <~gwpy.detector.channel.Channel>` for this
        `TabSummaryPlot`.

        :type: :class:`~gwpy.detector.channel.ChannelList`
        """
        return ChannelList(get_channel(c) for c in self._channels)

    @channels.setter
    def channels(self, channellist):
        self._channels = channellist

    @property
    def ifos(self):
        """Interferometer set for this `TabSummaryPlot`
        """
        return self.channels.ifos

    @property
    def tag(self):
        """File tag for this `DataPlot`.
        """
        try:
            return self._tag
        except AttributeError:
            state = re_cchar.sub('_', self.state is None and 'MULTI' or
                                      self.state.name)
            chans = "".join(map(str, self.channels))
            filts = "".join(map(str,
                [getattr(c, 'filter', getattr(c, 'frequency_response', ''))
                 for c in self.channels]))
            hash = hashlib.md5(chans+filts).hexdigest()[:6]
            type_ = re_cchar.sub('_', self.type)
            return '_'.join([state, hash, type_]).upper()

    @tag.setter
    def tag(self, filetag):
        self._tag = filetag or self.type.upper()

    @property
    def outputfile(self):
        """Output file for this `TabSummaryPlot`.
        """
        ifos = ''.join(self.ifos)
        tag = self.tag
        gps = floor(self.start)
        dur = ceil(self.end - self.start)
        return os.path.join(self._outdir,
                            '%s-%s-%d-%d.png' % (ifos, tag, gps, dur))

    @property
    def href(self):
        if self._href is None:
            return self.outputfile
        else:
            return self._href

    @href.setter
    def href(self, url):
        self._href = url and os.path.normpath(url) or None

    # ------------------------------------------------------------------------
    # TabSummaryPlot methods

    def parse_legend_kwargs(self, **defaults):
        """Pop the legend arguments from the `pargs` for this Plot
        """
        legendargs = defaults.copy()
        for key in self.pargs.keys():
            if re.match('legend[-_]', key):
                legendargs[key[7:]] = self.pargs.pop(key)
        return legendargs

    def parse_plot_kwargs(self, **defaults):
        """Pop keyword arguments for `Axes.plot` from the `pargs` for this Plot
        """
        plotargs = defaults.copy()
        plotargs.setdefault('label', self._parse_labels())
        chans_ = self.get_channel_groups()
        for kwarg in ['alpha', 'color', 'drawstyle', 'fillstyle', 'linestyle',
                      'linewidth', 'marker', 'markeredgecolor',
                      'markeredgewidth', 'markerfacecolor',
                      'markerfacecoloralt', 'markersize']:
            try:
                val = self.pargs.get(kwarg, self.pargs.get('%ss' % kwarg))
            except KeyError:
                continue
            else:
                plotargs[kwarg] = val
        chans = self.get_channel_groups().keys()
        for key, val in plotargs.iteritems():
            if (not isinstance(val, (list, tuple)) or len(val) != len(chans_)):
                plotargs[key] = [val]*len(self.get_channel_groups())
        out = []
        for i in range(len(chans)):
            out.append(dict((key, val[i]) for key, val in plotargs.items() if
                            val is not None and val[i] is not None))
        return out

    def _parse_labels(self, defaults=None):
        """Pop the labels for plotting from the `pargs` for this Plot
        """
        chans = self.get_channel_groups().keys()
        if defaults is None:
            defaults = chans
        labels = self.pargs.pop('labels', defaults)
        if isinstance(labels, (unicode, str)):
            labels = labels.split(',')
        labels = map(lambda s: rUNDERSCORE.sub(r'\_', str(s).strip('\n ')),
                     labels)
        while len(labels) < len(chans):
            labels.append(None)
        return labels

    def add_channel(self, channel):
        self._channels.append(channel)

    def get_channel_groups(self):
        """Find and group (mean, min, max) sets of channels
        for plotting.

        Returns
        -------
        groups : :class:`OrderedDict`
            dict of (channelname, channellist) pairs giving core channel
            name and an ordered list of channels. Ordering in preference
            of 'rms', 'mean', 'min', 'max'.
        """
        all_ = self.channels
        out = OrderedDict()
        for c in all_:
            name = c.texname.rsplit('.', 1)[0]
            if name in out.keys():
                out[name].append(c)
            else:
                out[name] = [c]
        order = ['rms', 'mean', 'min', 'max']
        for key in out.keys():
            out[key].sort(key=lambda c: c.name.split('.')[-1] in order and
                                        order.index(c.name.split('.')[-1])+1 or
                                        10)
        return out

    @classmethod
    def from_ini(cls, config, section, start, end, channels=None, **kwargs):
        """Define a new `DataPlot`.
        """
        # read parameters
        try:
            params = dict(config.nditems(section))
        except AttributeError:
            params = dict(config.items(section))
        for key in params:
            key2 = re_cchar.sub('_', key)
            if key != key2:
                params[key2] = params.pop(key)

        # get and check type
        ptype = re.sub('[\'\"]', '', params.pop('type'))
        if ptype != cls.type:
            warnings.warn("'%s' plot definition from configuration being "
                          "parsed by different plotting class '%s'"
                          % (ptype, cls.__name__))
        # get channels
        if channels is None:
            channels = params.pop('channels', [])
        if isinstance(channels, (unicode, str)):
            channels = split_channels(channels)
        # parse other parameters
        for key, val in params.iteritems():
            try:
                params[key] = eval(val)
            except NameError:
                pass
        params.update(kwargs)
        # format and return
        return cls(channels, start, end, **params)

    def queue(self, queue):
        """Submit this plot to the queue for processing.
        """
        # process data
        try:
            self.process()
        # announce done regardless of status
        finally:
            queue.get()
            queue.task_done()

    def process(self):
        """Process all data and generate the output file for this
        `SummaryPlot`.

        This function should be provided by all sub-classes, and should
        take no arguments.
        """
        raise NotImplementedError("This method should be provided by a "
                                  "sub-class")

    def finalize(self, outputfile=None):
        """Save the plot to disk and close.
        """
        # quick fix for x-axis labels hitting the axis
        for ax in self.plot.axes:
            ax.tick_params(axis='x', pad=10)
            ax.xaxis.labelpad = 10
        # save figure and close
        if outputfile is None:
            outputfile = self.outputfile
        self.plot.save(outputfile)
        self.plot.close()
        return outputfile

register_plot(DataPlot)
