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

import os.path
import re
import warnings
from collections import OrderedDict
from math import (floor, ceil)
from numbers import Number

from six import string_types
from six.moves.urllib.parse import urlparse

from matplotlib import (rcParams, rc_context)

from gwpy.segments import Segment
from gwpy.detector import ChannelList
from gwpy.plot import Plot

from ..channels import (get_channel, split as split_channels,
                        split_combination as split_channel_combination)
from ..config import GWSummConfigParser
from ..state import get_state
from ..utils import (vprint, safe_eval, re_quote)
from . import utils as putils
from .registry import register_plot

__all__ = ['SummaryPlot', 'DataPlot']

re_cchar = re.compile(r"[\W\s_]+")

putils.AXES_PARAMS.extend([
    'insetlabels',  # for segment plotting
])
NON_PLOT_PARAMS = set(putils.FIGURE_PARAMS + putils.AXES_PARAMS)


# -- utilities ----------------------------------------------------------------

def format_label(label):
    label = str(label).strip('\n ')
    label = re_quote.sub('', label)
    return putils.usetex_tex(label)


# -- basic Plot object --------------------------------------------------------

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

    def __init__(self, href=None, src=None, new=True, caption=''):
        self.href = href
        if src:
            self.src = src
        self.new = new
        self.caption = caption

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

    @property
    def caption(self):
        """HTML <fancybox plot> title attribute for this `SummaryPlot`.
        """
        return self._caption

    @caption.setter
    def caption(self, text):
        self._caption = text

    # -- SummaryPlot methods -----------------------

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
        a descriptive tag for this `DataPlot`, used as part of the output
        file name
    outdir : `str`
        output directory path for this `DataPlot`, defaults to
        the current directory
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
    #: name for DataPlot subclass
    type = 'data'
    #: plot() call style
    _single_call = False
    #: dict of default plotting kwargs
    defaults = {}
    #: list of parameters parsed for `plot()` calls
    DRAW_PARAMS = list(putils.ARTIST_PARAMS)

    def __init__(self, channels, start, end, state=None, outdir='.',
                 tag=None, pid=None, href=None, new=True, all_data=False,
                 read=True, fileformat='png', caption=None, **pargs):
        super(DataPlot, self).__init__(href=href, new=new, caption=caption)
        if isinstance(channels, string_types):
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

    # -- properties -----------------------------

    @property
    def span(self):
        """The GPS [start, stop) interval for this `DataPlot`.
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
        if isinstance(state_, string_types):
            self._state = get_state(state_)
        else:
            self._state = state_

    @property
    def channels(self):
        """List of data-source
        :class:`Channels <~gwpy.detector.channel.Channel>` for this
        `DataPlot`.

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
        return type(self.channels)(OrderedDict.fromkeys(
            c2 for c in self.channels for
            c2 in split_channel_combination(c.ndsname)
        ).keys())

    @property
    def ifos(self):
        """Interferometer set for this `DataPlot`
        """
        return set([c.ifo for c in self.allchannels if c.ifo])

    @property
    def tag(self):
        """File tag for this `DataPlot`.
        """
        try:
            return self._tag
        except AttributeError:
            state = re_cchar.sub(
                '_',
                self.state is None and 'MULTI' or self.state.name).rstrip('_')
            type_ = re_cchar.sub('_', self.type)
            self._tag = '_'.join([state, self.pid, type_]).upper()
            return self.tag

    @tag.setter
    def tag(self, filetag):
        if filetag is None:
            del self.tag
        else:
            self._tag = filetag

    @tag.deleter
    def tag(self):
        del self._tag

    @property
    def pid(self):
        try:
            return self._pid
        except AttributeError:
            chans = "".join(map(str, self.channels))
            filts = "".join(map(str, [
                getattr(c, 'filter', getattr(c, 'frequency_response', ''))
                for c in self.channels]))
            self._pid = putils.hash(chans + filts)
            return self.pid

    @pid.setter
    def pid(self, id_):
        self._pid = str(id_)

    @pid.deleter
    def pid(self):
        del self._pid

    @property
    def outputfile(self):
        """Output file for this `DataPlot`.
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

    # -- read-only plot properties --------------

    def _is_log(self, axis):
        scale = '{0}scale'.format(axis)
        try:
            return self.pargs[scale] == 'log'
        except KeyError:
            try:
                return self.pargs['log{0}'.format(axis)]
            except KeyError:
                try:
                    ax = self.plot.gca()
                except AttributeError:  # plot not generated yet
                    return False
                return getattr(ax, 'get_{0}'.format(scale))() == 'log'

    @property
    def logx(self):
        return self._is_log('x')

    @property
    def logy(self):
        return self._is_log('y')

    # -- basic methods --------------------------

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
            if c.ifo == 'G1' and re.search(r'-(av|min|max)\Z', c.texname):
                name = c.texname.rsplit('-', 1)[0]
            else:
                name = c.texname.rsplit('.', 1)[0]
            if ' ' in c.texname:
                out.append((c.texname, [c]))
            else:
                try:
                    id_ = list(zip(*out))[0].index(name)
                except (IndexError, ValueError):
                    out.append((name, [c]))
                else:
                    out[id_][1].append(c)
        order = ['rms', 'mean', 'av', 'min', 'max']
        for channel, clist in out:
            clist.sort(key=lambda c: c.name.split('.')[-1] in order and
                       order.index(c.name.split('.')[-1])+1 or 10)
        return out

    @classmethod
    def from_ini(cls, config, section, start, end, channels=None, **kwargs):
        """Define a new `DataPlot`.
        """
        config = GWSummConfigParser.from_configparser(config)

        # read parameters
        try:
            params = dict(config.nditems(section))
        except AttributeError:
            params = dict(config.items(section))

        # get and check type
        ptype = re.sub(r'[\'\"]', '', params.pop('type'))
        if ptype != cls.type:
            warnings.warn("'%s' plot definition from configuration being "
                          "parsed by different plotting class '%s'"
                          % (ptype, cls.__name__))
        # get channels
        if channels is None:
            channels = params.pop('channels', [])
        if isinstance(channels, string_types):
            channels = split_channels(channels)
        # parse specific parameters
        if 'all-data' in params:
            params['all_data'] = params.pop('all-data')

        # parse other parameters
        for key, val in params.items():
            params[key] = safe_eval(val)
        params.update(kwargs)

        # format and return
        return cls(channels, start, end, **params)

    # -- plot parameter parsing -----------------

    def _parse_param(self, pdict, key, allow_plural=False):
        """Parse a configuration parameter for this Plot from a dict

        Parameters
        ----------
        pdict : `dict`
            the dict to evaluate from

        key : `str`
            the key to evaluated

        allow_plural : `bool`, optional
            try to find plural version of ``key`` if singular not found,
            i.e. ``'colors'`` instead of ``'color'``, default is `False`

        Returns
        -------
        values : `list`
            a mapping of the parsed value for each channel in this plot

        Raises
        ------
        KeyError
            if ``key`` is not found in ``pdict``

        Notes
        -----
        This method uses `gwsumm.utils.safe_eval` to `eval` strings into
        python objects.

        Examples
        --------
        Consider a `DataPlot` with two channels to display

        >>> from gwsumm.plot import DataPlot
        >>> a = DataPlot(['channel1', 'channel2'], 0, 1)

        We can then parse params as follows:

        >>> a._parse_param({'marker': 'x'}, 'marker')
        ['x', 'x']
        >>> a._parse_param({'linestyles': ['-', '--']}, 'linestyle'}
        ['-', '--']
        >>> a._parse_param({'colors': "'red','green'"}, 'color'}
        ['red', 'green']
        """
        # parse keyword
        try:
            val = pdict.pop(key)
        except KeyError as e:
            if not allow_plural:
                raise
            # check for plural
            try:
                val = pdict.pop('%ss' % key)
            except KeyError:
                raise e

        # evaluate (safely) allowing references to self as 'plot'
        if isinstance(val, string_types):
            try:
                val = safe_eval(val, locals_={'plot': self, 'self': self})
            except ZeroDivisionError:  # e.g. zero livetime
                val = 0

        # don't use sets
        if isinstance(val, set):
            val = list(val)

        # if a single-call style plot, just return the value as given
        if self._single_call:
            return val

        # otherwise convert to a 1<->1 mapping with the channels list
        nchans = len(self.get_channel_groups())
        if not isinstance(val, (list, tuple)) or len(val) != nchans:
            return [val] * nchans

        return val

    def _parse_extra_params(self, prefix, **defaults):
        """Parse parameters for an extra plot element

        Parameters
        ----------
        prefix : `str`
            the text prefix identifying parameters for the extra element

        **defaults
            any default options to use

        Returns
        -------
        params : `dict`
        """
        re_prefix = re.compile(r'\A%s[-_]' % prefix.rstrip('-_'))
        extras = defaults.copy()
        for key in list(self.pargs):
            m = re_prefix.match(key)
            if m:
                extras[key[m.span()[1]:]] = safe_eval(self.pargs.pop(key))
        return extras

    def parse_legend_kwargs(self, **defaults):
        """Pop the legend arguments from the `pargs` for this `Plot`
        """
        return self._parse_extra_params('legend', **defaults)

    def parse_plot_kwargs(self, **defaults):
        """Pop keyword arguments for `Axes.plot` from the `pargs` for this Plot
        """
        plotargs = defaults.copy()
        plotargs.setdefault('label', self._parse_labels())

        # loop over known Axes.plot kwargs and parse
        for kwarg in self.DRAW_PARAMS:
            try:
                plotargs[kwarg] = self._parse_param(
                    self.pargs, kwarg, allow_plural=True)
            except KeyError:
                pass

        # normalise log scale parameters
        # TODO: this can be removed once existing config files have been
        #       updated to not use logx and logy options
        for axis in ('x', 'y'):
            logp = 'log{}'.format(axis)
            try:
                log = self.pargs.pop(logp)
            except KeyError:
                continue
            if log is not self._is_log(axis):
                scale = 'log' if log else 'linear'
                self.pargs['{}scale'.format(axis)] = scale

        # if this plot is a single-call plot (where all objects get plotted
        # in a single call out to a ax.plot()-style method) just return
        # the params as a dict of lists
        if self._single_call:
            return plotargs

        # otherwise, map to a list of dicts (one per channel)
        out = []
        nchans = len(self.get_channel_groups())
        for i in range(nchans):
            out.append(dict((key, val[i]) for key, val in plotargs.items() if
                            val is not None and val[i] is not None))
        return out

    def _parse_labels(self, defaults=None):
        """Pop the labels for plotting from the `pargs` for this Plot
        """
        # set default label to show channel name
        chans = list(zip(*self.get_channel_groups()))[0]
        if defaults is None:
            defaults = chans

        # parse user labels
        try:
            labels = self._parse_param(self.pargs, 'label', allow_plural=True)
        except KeyError:
            labels = defaults

        if list(set(labels)) == [labels[0]] and labels[0] is not None:
            labels = labels[0].split(',')

        # escape underscores
        labels = list(map(format_label, labels))

        # fill gaps with None
        while len(labels) < len(chans):
            labels.append(None)

        return labels

    def parse_rcParams(self, params):
        """Parse matplotlib rcParams settings from a dict of plot params
        """
        self.rcParams = {}
        for key in list(params):
            if key in rcParams:
                self.rcParams[key] = safe_eval(params.pop(key))
        return self.rcParams

    def parse_list(self, prefix, **defaults):
        """Parse a list of something from parameters

        This enables listing `hline`s (for example) in the config as

        [plot-blah]
        hline = 100
        hline-linestyle = '--'
        hline-color = 'red'
        hline2 = 200
        hline2-linestyle = '--'
        hline2-color = 'blue'

        Returns an `OrderedDict` with keys matching the primary parsed
        value, and values as everything else, e.g.

        {100: {'linestyle': '--', 'color': 'red'},
         200: {'linestyle': '--', 'color': 'blue'},}
        """
        items = OrderedDict()
        re_prefix = re.compile(r'{0}(\d+)?\Z'.format(prefix))
        keys = sorted(list(self.pargs))
        while True:
            for i, key in enumerate(keys):
                if re_prefix.match(key):
                    primary = safe_eval(self.pargs.pop(key))
                    items[primary] = self._parse_extra_params(key, **defaults)
                    break  # go to next iteration
            else:  # no references matched, so stop
                break
            keys = sorted(list(self.pargs)[i:])

        return items

    # -- figure processing ----------------------

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

    def init_plot(self, data=[], FigureClass=Plot, geometry=(1, 1),
                  projection='rectilinear', sharex=True, sharey=True,
                  **kwargs):
        """Initialise the Figure and Axes objects for this `DataPlot`.
        """
        # update plot defaults using channel data
        self._update_defaults_from_channels()

        # strip figure and axes params from pargs
        for key in NON_PLOT_PARAMS:
            try:
                kwargs.setdefault(key, self.pargs.pop(key))
            except KeyError:
                continue
            # escape text for TeX
            if key in ('title', 'xlabel', 'ylabel'):
                kwargs[key] = putils.usetex_tex(kwargs[key])

        # create figure
        self.plot = FigureClass(*data, geometry=geometry,
                                projection=projection, sharex=sharex,
                                sharey=sharey, **kwargs)
        return self.plot

    def _update_defaults_from_channels(self):
        """Update default plotting params from channel attributes

        This method is called at the start of DataPlot.init_plot(),
        so should be populated in any subclasses that want it
        """
        pass

    def finalize(self, outputfile=None, close=True, **savekwargs):
        """Save the plot to disk and close.
        """
        # save figure and close (build both png and pdf for pdf choice)
        if outputfile is None:
            outputfile = self.outputfile
        if not isinstance(outputfile, string_types):
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
            except (IOError, RuntimeError, IndexError) as e:
                warnings.warn("Caught %s: %s [retrying...]"
                              % (type(e).__name__, str(e)))
                self.plot.save(fp, **savekwargs)
            if isinstance(fp, string_types):
                vprint("        %s written\n" % fp)
        if close:
            self.plot.close()
        return outputfile

    def apply_parameters(self, *axes, **pargs):
        keys = sorted(list(pargs),
                      key=lambda x: 1 if x in ('xscale', 'yscale') else 2)
        for ax in axes:
            for key in keys:
                if key.startswith('no-'):  # skip no-xxx keys
                    continue
                val = pargs[key]
                if key in ['xlim', 'ylim'] and isinstance(val, string_types):
                    val = eval(val)
                elif key == 'grid':
                    self._apply_grid_params(ax, val)
                    continue
                try:
                    getattr(ax, 'set_%s' % key)(val)
                except AttributeError:
                    setattr(ax, key, val)

    def _apply_grid_params(self, ax, val):
        if val in ('major', 'minor'):
            ax.grid(True, which=val)
        else:
            ax.grid(val)

    def add_hvlines(self):
        """Add horizontal and vertical lines to this `DataPlot`

        These should be defined in the configuration via the `hline` and
        `vline` keys.
        """
        for key in ('hline', 'vline'):
            lines = self.parse_list(key, linestyle='--', color='red')
            for ax in self.plot.axes:
                axline = getattr(ax, 'ax{0}'.format(key))
                for val, params in lines.items():
                    if isinstance(val, Number):
                        val = [val]
                    for x in val:
                        axline(x, **params)


register_plot(DataPlot)


# -- custom plot types --------------------------------------------------------

class BarPlot(DataPlot):
    """`DataPlot` with bars
    """
    type = 'bar'
    DRAW_PARAMS = ['width', 'bottom', 'color', 'edgecolor', 'linewidth',
                   'xerr', 'yerr', 'ecolor', 'capsize', 'error_kw',
                   'align', 'orientation', 'log', 'alpha', 'rasterized']


class PiePlot(DataPlot):
    type = 'pie'
    DRAW_PARAMS = ['explode', 'colors', 'autopct', 'pctdistance', 'shadow',
                   'labeldistance', 'startangle', 'radius', 'counterclock',
                   'wedgeprops', 'textprops', 'rasterized']
