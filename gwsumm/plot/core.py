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

import os
import re
from math import (floor, ceil)

try:
    from collections import OrderedDict
except ImportError:
    from astropy.utils import OrderedDict

from gwpy.detector import (Channel, ChannelList)

from .. import globalv

__all__ = ['TabPlot', 'PlotList']

re_cchar = re.compile("[\W\s_]+")


class TabPlot(object):
    """Defines a plot to be displayed on a `Tab`

    This object is designed to be sub-classed by specific plot types
    to be registered to a given name, e.g. the TimeSeriesTabPlot class
    is automatically registered to the 'timeseries' name.

    Parameters
    ----------
    channels : `list`
        a list of channel names that define the data sources for this
        `TabPlot`
    segment : :class:`~gwpy.segments.segments.Segment`
        the defining GPS [start, stop) `Segment` for this `TabPlot`.
         This is used to constructed the LIGO-T050017 compliant output
         file name for this plot
    tag : `str`
        a descriptive tag for this `TabPlot`, used as part of the output
        file name
    outdir : `str`
        output directory path for this `TabPlot`, defaults to the current
        directory
    href : `str`
        custom URL for this plot to link towards.
    **kwargs
        all other keyword arguments to be passed to this plot's
        :meth:`process` method.

    Attributes
    ----------
    channels
    state
    gpsstart
    gpsend
    outputfile
    href

    Notes
    -----
    All sub-classes of this object must provide the following methods

    ==============  ===================================================
    :math`process`  to method by which the data are loaded/generated \
                    and the plot actually written
    ==============  ===================================================
    """
    #: Figure subclass
    FigureClass = None
    #: Axes subclass
    AxesClass = None
    #: name for TabPlot subclass
    type = None
    #: dict of default plotting kwargs
    defaults = {}

    def __init__(self, channels, state='all', outdir='.', href=None, **kwargs):
        self.channels = channels
        self.state = state
        self.outdir = outdir
        self.plotargs = kwargs
        self.href = href

    # ------------------------------------------------------------------------
    # TabPlot properties

    @property
    def channels(self):
        """List of data-source
        :class:`Channels <~gwpy.detector.channel.Channel>` for this
        `TabPlot`.

        :type: :class:`~gwpy.detector.channel.ChannelList`
        """
        return self._channels

    @channels.setter
    def channels(self, channellist):
        self._channels = ChannelList()
        for ch in channellist:
            if not isinstance(ch, Channel):
                ch = Channel(ch)
            self._channels.append(ch)

    @property
    def state(self):
        """Instrumental :class:`~gwsumm.state.SummaryState` defining
        the useful data segments for this `TabPlot`
        """
        return self._state

    @state.setter
    def state(self, state_):
        if isinstance(state_, basestring):
            self._state = globalv.STATES[state_]
        else:
            self._state = state_

    @property
    def gpsstart(self):
        """GPS start time for this `TabPlot`
        """
        return self.state.extent[0]

    @property
    def gpsend(self):
        """GPS end time for this `TabPlot`
        """
        return self.state.extent[1]

    @property
    def href(self):
        """HTML <a> URL for this `TabPlot`.
        """
        if self._href:
            return self._href
        else:
            return self.outputfile

    @href.setter
    def href(self, url):
        self._href = url

    @property
    def tag(self):
        """Descriptive filetag for this `TabPlot`.
        """
        raise NotImplementedError("Sub-classes must define a tag")

    @property
    def ifos(self):
        """Interferometer set for this `TabPlot`
        """
        return self.channels.ifos

    @property
    def outputfile(self):
        """Output file path for this `TabPlot`.
        """
        ifos = ''.join(self.ifos)
        desc = '_'.join([re_cchar.sub('_', self.state.name.upper()),
                         self.tag.upper(),
                         self.type.upper()])
        return os.path.join(self.outdir, '%s-%s-%d-%d.png'
                                         % (ifos, desc, floor(self.gpsstart),
                                            ceil(abs(self.state.extent))))

    # ------------------------------------------------------------------------
    # TabPlot methods

    def process(self):
        """Generate a plot using the data pre-loaded into the
        globalv memory cache

        This method should be sub-classed by other plot types
        """
        raise NotImplementedError("Actually generating plots hasn't been "
                                  "implemented yet")

    def find_mmm_channels(self):
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
        all_.sort(key=lambda c: str(c))
        for c in all_:
            name = str(c).split('.')[0]
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

    # ------------------------------------------------------------------------
    # TabPlot comparisons

    def __eq__(self, other):
        """Compare this `TabPlot` to another.

        Returns
        -------
        `True`
            if the metadata for ``self`` match those of ``other``
        `False`
            otherwise
        """
        if not isinstance(other, self.__class__):
            return False
        if not self.channels == other.channels:
            return False
        if not self.type == other.type:
            return False
        if not self.outdir == other.outdir:
            return False
        return True

    def __repr__(self):
        return '<%sPlotRequest(%s, %s)>' % (self.type.title(), self.channels,
                                            self.outdir)


class PlotList(list):
    """A searchable list of plots that have been configured to be
    generated
    """

    @property
    def timeseries(self):
        return self.find_all('timeseries')

    @property
    def spectra(self):
        return self.find_all('spectrum')

    @property
    def spectrograms(self):
        return self.find_all('spectrogram')

    @property
    def segments(self):
        return self.find_all('segments')

    def find_all(self, type):
        return [s for s in self if s.type.lower() == type]

    def unique(self):
        new = self.__class__()
        for entry in self:
            if not entry in new:
                new.append(entry)
        return new

