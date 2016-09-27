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

from matplotlib import rc_context

from gwpy.segments import Segment
from gwpy.detector import (Channel, ChannelList)
from gwpy.plotter.utils import rUNDERSCORE

from . import rcParams
from .registry import register_plot
from .. import globalv
from ..channels import (get_channel, split as split_channels)
from ..utils import (vprint, safe_eval)

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

    def __init__(self, href=None, src=None, new=True):
        self.href = href
        if src:
            self.src = src
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

    @property
    def src(self):
        try:
            return self._src
        except AttributeError:
            return self.href

    @src.setter
    def src(self, url):
        self._src = url

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
    #: list of parameters parsed for `plot()` calls
    DRAW_PARAMS = ['alpha', 'color', 'drawstyle', 'fillstyle', 'linestyle',
                   'linewidth', 'marker', 'markeredgecolor',
                   'markeredgewidth', 'markerfacecolor',
                   'markerfacecoloralt', 'markersize', 'valid', 'edgecolor',
                   'bins', 'range', 'normed', 'weights', 'cumulative',
                   'bottom', 'histtype', 'align', 'orientation', 'rwidth',
                   'log', 'stacked', 'logbins', 'linecolor',
                   'facecolor', 'rasterized']

    def __init__(self, channels, start, end, state=None, outdir='.',
                 tag=None, pid=None, href=None, new=True, all_data=False,
                 read=True, fileformat='png', **pargs):
        super(DataPlot, self).__init__(href=href, new=new)
        if isinstance(channels, str):
            channels = split_channels(channels)
        self.channels = channels
        self.span = (start, end)
        self.state = state
        self._outdir = outdir
        if tag is not None:
            self.tag = tag
        if pid is not None:
            self.pid = pid
        # allow user to specify all-data instead of all_data as kwarg
        # mainly for INI-parsing convenience, should fix properly
        self.all_data = pargs.pop('all-data', all_data)
        self.pargs = self.defaults.copy()
        self.pargs.update(pargs)
        self.parse_rcParams(self.pargs)
        self.plot = None
        self.read = read
        self.fileformat = fileformat

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
    def allchannels(self):
        """List of all unique channels for this plot
        """
        out = type(self.channels)()
        for c in self.channels:
            for m in Channel.MATCH.finditer(c.ndsname):
                c2 = get_channel(c.ndsname[m.start():m.end()])
                if c2 not in out:
                    out.append(c2)
        return out

    @property
    def ifos(self):
        """Interferometer set for this `TabSummaryPlot`
        """
        return set([c.ifo for c in self.allchannels if c.ifo])

    @property
    def tag(self):
        """File tag for this `DataPlot`.
        """
        try:
            return self._tag
        except AttributeError:
            state = re_cchar.sub('_', self.state is None and 'MULTI' or
                                      self.state.name).rstrip('_')
            type_ = re_cchar.sub('_', self.type)
            self._tag = '_'.join([state, self.pid, type_]).upper()
            return self.tag

    @tag.setter
    def tag(self, filetag):
        self._tag = filetag or self.type.upper()

    @tag.deleter
    def tag(self):
        del self._tag

    @property
    def pid(self):
        try:
            return self._pid
        except:
            chans = "".join(map(str, self.channels))
            filts = "".join(map(str,
                [getattr(c, 'filter', getattr(c, 'frequency_response', ''))
                 for c in self.channels]))
            self._pid = hashlib.md5(chans+filts).hexdigest()[:6]
            return self.pid

    @pid.setter
    def pid(self, id_):
        self._pid = str(id_)

    @pid.deleter
    def pid(self):
        del self._pid

    @property
    def outputfile(self):
        """Output file for this `TabSummaryPlot`.
        """
        ifos = ''.join(sorted(self.ifos))
        tag = self.tag
        gps = floor(self.start)
        dur = ceil(self.end - self.start)
        return os.path.join(
            self._outdir,
            '%s-%s-%d-%d.%s' % (ifos, tag, gps, dur, self.fileformat))

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
        for kwarg in self.DRAW_PARAMS:
            try:
                val = self.pargs.pop(kwarg)
            except KeyError:
                try:
                    val = self.pargs.pop('%ss' % kwarg)
                except KeyError:
                    val = None
            if val is not None and isinstance(val, str) and 'self' in val:
                try:
                    plotargs[kwarg] = safe_eval(val, locals_={'self': self})
                except ZeroDivisionError:
                    plotargs[kwarg] = 0
            elif val is not None:
                plotargs[kwarg] = safe_eval(val)
        chans = zip(*self.get_channel_groups())[0]
        for key, val in plotargs.iteritems():
            if (key.endswith('color') and isinstance(val, (list, tuple)) and
                    isinstance(val[0], (int, float))):
                plotargs[key] = [val]*len(self.get_channel_groups())
            elif (not isinstance(val, (list, tuple)) or len(val) != len(chans)):
                plotargs[key] = [val]*len(self.get_channel_groups())
        out = []
        for i in range(len(chans)):
            out.append(dict((key, val[i]) for key, val in plotargs.items() if
                            val is not None and val[i] is not None))
        return out

    def _parse_labels(self, defaults=None):
        """Pop the labels for plotting from the `pargs` for this Plot
        """
        chans = zip(*self.get_channel_groups())[0]
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

    def parse_rcParams(self, params):
        """Parse matplotlib rcParams settings from a dict of plot params
        """
        self.rcParams = {}
        for key in params.keys():
            if key in rcParams:
                self.rcParams[key] = params.pop(key)
        return self.rcParams

    def add_channel(self, channel):
        self._channels.append(channel)

    def get_channel_groups(self):
        """Find and group (mean, min, max) sets of channels for plotting.

        Returns
        -------
        groups : `list` of `tuple`
            list of (channelname, channellist) tuples giving core channel
            name and an ordered list of channels. Ordering in preference
            of 'rms', 'mean', 'min', 'max'.

        Notes
        -----
        This method used to return an `OrderedDict`, but was changed to
        return a `list` of `tuple` to enable plotting a channel multiple
        times on a plot, for whatever reason.
        """
        all_ = self.channels
        out = []
        for c in all_:
            if c.ifo == 'G1' and re.search('-(av|min|max)\Z', c.texname):
                name = c.texname.rsplit('-', 1)[0]
            else:
                name = c.texname.rsplit('.', 1)[0]
            if ' ' in c.texname:
                out.append((c.texname, [c]))
            else:
                try:
                    id_ = zip(*out)[0].index(name)
                except (IndexError, ValueError):
                    out.append((name, [c]))
                else:
                    out[id_][1].append(c)
        order = ['rms', 'mean', 'av', 'min', 'max']
        for channel, clist in out:
            clist.sort(key=lambda c: c.name.split('.')[-1] in order and
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
        # parse specific parameters
        if 'all-data' in params:
            params['all_data'] = params.pop('all-data')

        # parse other parameters
        for key, val in params.iteritems():
            params[key] = safe_eval(val)
        params.update(kwargs)
        # escape text
        for key, val in params.iteritems():
            if key in ['title', 'ylabel', 'xlabel']:
                params[key] = re.sub(r'(?<!\\)_(?!.*{)', '\_', params[key])
        # format and return
        return cls(channels, start, end, **params)

    # -- figure processing ----------------------------------------------------

    def process(self, outputfile=None, close=True):
        with rc_context(rc=self.rcParams):
            return self.draw()

    def draw(self):
        """Process all data and generate the output file for this
        `SummaryPlot`.

        This function should be provided by all sub-classes, and should
        take no arguments.
        """
        raise NotImplementedError("This method should be provided by a "
                                  "sub-class")

    def finalize(self, outputfile=None, close=True, **savekwargs):
        """Save the plot to disk and close.
        """
        # customise axes
        for ax in self.plot.axes:
            # quick fix for x-axis labels hitting the axis
            if not self.type == 'bar' or self.type.endswith('-bar'):
                ax.tick_params(axis='x', pad=10)
                ax.xaxis.labelpad = 10
            # move title up to create gap between axes
            if ax.get_title() and ax.title.get_position()[1] == 1.0:
                ax.title.set_y(1.01)
            # lighten color of axes and legend borders
            color = rcParams['grid.color']
            for edge in ax.spines:
                ax.spines[edge].set_edgecolor(color)
            if ax.legend_ and ax.legend_.get_frame().get_edgecolor() != 'none':
                ax.legend_.get_frame().set_edgecolor(color)
        # customise colorbars
        for cb in self.plot.colorbars:
            cb.outline.set_edgecolor(color)
        # save figure and close (build both png and pdf for pdf choice)
        if outputfile is None:
            outputfile = self.outputfile
        if not isinstance(outputfile, str):
            extensions = [None]
        elif outputfile.endswith('.pdf'):
            extensions = ['.pdf', '.png']
        else:
            extensions = [os.path.splitext(outputfile)[1]]

        for ext in extensions:
            try:
                fp = '%s%s' % (os.path.splitext(outputfile)[0], ext)
            except AttributeError:
                fp = outputfile
            try:
                self.plot.save(fp, **savekwargs)
            except (IOError, RuntimeError) as e:
                warnings.warn("Caught %s: %s [retrying...]"
                             % (type(e).__name__, str(e)))
                self.plot.save(fp, **savekwargs)
            if isinstance(fp, str):
                vprint("        %s written\n" % fp)
        if close:
            self.plot.close()
        return outputfile

    def apply_parameters(self, ax, **pargs):
        for key in pargs:
            val = pargs[key]
            if key in ['xlim', 'ylim'] and isinstance(val, str):
                val = eval(val)
            try:
                getattr(ax, 'set_%s' % key)(val)
            except AttributeError:
                setattr(ax, key, val)


register_plot(DataPlot)


# -- custom plot types --------------------------------------------------------


class _SingleCallPlot(object):
    """Custom plot mixin to parse plot kwargs for a single call

    """
    DRAW_PARAMS = []

    def parse_plot_kwargs(self, defaults=dict()):
        """Parse pie() keyword arguments
        """
        plotargs = defaults.copy()
        plotargs.setdefault('labels', self._parse_labels())
        for kwarg in self.DRAW_PARAMS:
            try:
                val = self.pargs.pop(kwarg)
            except KeyError:
                try:
                    val = self.pargs.pop('%ss' % kwarg)
                except KeyError:
                    val = None
            if val is not None:
                plotargs[kwarg] = safe_eval(val)
                plotargs[kwarg] = val
        return plotargs



class BarPlot(_SingleCallPlot, DataPlot):
    """`DataPlot` with bars
    """
    type = 'bar'
    DRAW_PARAMS = ['width', 'bottom', 'color', 'edgecolor', 'linewidth',
                   'xerr', 'yerr', 'ecolor', 'capsize', 'error_kw',
                   'align', 'orientation', 'log', 'alpha', 'rasterized']


class PiePlot(_SingleCallPlot, DataPlot):
    type = 'pie'
    DRAW_PARAMS = ['explode', 'colors', 'autopct', 'pctdistance', 'shadow',
                   'labeldistance', 'startangle', 'radius', 'counterclock',
                   'wedgeprops', 'textprops', 'rasterized']
