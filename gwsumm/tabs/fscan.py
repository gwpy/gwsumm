# coding=utf-8
# Copyright (C) Duncan Macleod (2013)
#
# This file is part of GWSumm
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
# along with GWSumm.  If not, see <http://www.gnu.org/licenses/>

"""Custom `SummaryTab` for the output of the FScan algorithm.
"""

import os
import glob
from dateutil import parser

from MarkupPy import markup

from .registry import (get_tab, register_tab)

from gwdetchar.io import html

from .. plot import get_plot
from ..mode import Mode
from ..config import GWSummConfigParser

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__all__ = ['FscanTab']

base = get_tab('default')
SummaryPlot = get_plot(None)


class FscanTab(base):
    """Custom tab displaying a summary of Fscan results.
    """
    type = 'fscan'

    def __init__(self, *args, **kwargs):
        if kwargs['mode'] != Mode.day:
            raise RuntimeError("FscanTab is only available in %s mode."
                               % Mode.day.name)
        super(FscanTab, self).__init__(*args, **kwargs)

    @classmethod
    def from_ini(cls, config, section, **kwargs):
        """Define a new `FscanTab` from a `ConfigParser`.
        """
        # parse generic configuration
        new = super(FscanTab, cls).from_ini(config, section, **kwargs)
        new.set_layout([2])

        # work out day directory and url
        new.channel = config.get(section, 'channel')
        base = os.path.normpath(config.get(section, 'base-directory'))
        dirs = glob.glob(base)
        if len(dirs) == 0:
            d = None
        elif len(dirs) > 1:
            d = dirs
        else:
            d = dirs[0]
        if d is not None:
            new.directory = os.path.join(d, new.channel.replace(':', '_'))
        else:
            new.directory = None

        # get navigation urls
        new.navigation = []
        for key, val in filter(lambda x: x[0].startswith('navigation'),
                               config.items(section)):
            system = key.split('-', 1)[1]
            new.navigation.append((system, val))

        return new

    def process(self, config=GWSummConfigParser(), **kwargs):
        # find all plots
        self.plots = []
        self.line_count_plots = []
        if isinstance(self.directory, str):
            plots = sorted(
                glob.glob(os.path.join(self.directory, 'spec_*.png')),
                key=lambda p: float(os.path.basename(p).split('_')[1]))
            for p in plots:
                home_, postbase = p.split('/public_html/', 1)
                user = os.path.split(home_)[1]
                if 'line_count' not in p:
                    self.plots.append(SummaryPlot(
                        href='/~%s/%s' % (user, postbase)))
                else:
                    self.line_count_plots.append(SummaryPlot(
                        href='/~%s/%s' % (user, postbase)))
                if ('line_count' not in p and
                        '0.00_100.00' in p and
                        '_2.png' in p):
                    self.line_count_plots.append(SummaryPlot(
                        href='/~%s/%s' % (user, postbase)))

    def write_state_html(self, state):
        """Write the '#main' HTML content for this `FscanTab`.
        """
        page = markup.page()

        if self.directory is None:
            page.add(html.alert((
                "No analysis was performed for this period, "
                "please try again later.",
                "If you believe these data should have been found, please "
                "contact %s."
                % markup.oneliner.a('the DetChar group',
                                    class_='alert-link',
                                    href='mailto:detchar@ligo.org'),
            ), context='warning', dismiss=False))
            self.directory = []
        elif isinstance(self.directory, list):
            page.add(html.alert((
                "Multiple results directories were found for this period, "
                "cannot determine correct directory.",
                "If you believe this to be an error, please contact %s."
                % markup.oneliner.a('the DetChar group',
                                    class_='alert-link',
                                    href='mailto:detchar@ligo.org'),
            ), context='warning', dismiss=False))
            page.hr(class_='row-divider')
        elif not self.plots:
            page.add(html.alert((
                "This analysis produced no plots.",
                "If you believe these data should have been found, please "
                "contact %s."
                % markup.oneliner.a('the DetChar group',
                                    class_='alert-link',
                                    href='mailto:detchar@ligo.org'),
            ), context='warning', dismiss=False))
            page.hr(class_='row-divider')
            self.directory = [self.directory]
        else:
            self.directory = [self.directory]

        # link full results
        for d in self.directory:
            # parse date
            datestr = d.split(os.path.sep)[-2].split('_')[1:]
            datestr = '{0}-{1}-{2} {3}:{4}:{5} {6}'.format(*datestr)
            date = parser.parse(datestr).strftime('%Y-%m-%d %H:%M:%S %Z')
            # find HTML
            try:
                index = glob.glob(os.path.join(
                    d, self.channel.replace(':', '_', 1), '*html'))[0]
            except IndexError:
                index = d + os.path.sep
            home_, postbase = index.split('/public_html/', 1)
            user = os.path.split(home_)[1]
            index = '/~%s/%s' % (user, postbase)
            page.div(class_='btn-group')
            page.a('Full Fscan results for %s' % date,
                   href=index, rel='external', target='_blank',
                   class_='btn btn-info btn-xl')
            page.div.close()

        if self.line_count_plots:
            page.h2('', class_='mt-4 mb-2')
            page.div(class_='card border-light card-body scaffold shadow-sm')
            plt1 = self.line_count_plots[::3]
            plt2 = self.line_count_plots[1::3]
            plt3 = self.line_count_plots[2::3]
            page.p('Line count figure of merit (0 - 100 Hz):')
            for triple in list(zip(plt1, plt2, plt3)):
                page.div(class_='row')
                for p in triple:
                    page.div(class_="col-sm-4 mb-10")
                    page.a(href=p.href, class_="fancybox plot",
                           **{'data-fancybox-group': 1})
                    page.img(class_='img-fluid w-100', src=p.href)
                    page.a.close()
                    page.div.close()
                page.div.close()
            page.div.close()

        if self.navigation:
            page.h2('Fscan links:', class_='mt-4 mb-2')
        for i, (key, url) in enumerate(self.navigation):
            if not i % 4:
                page.div(class_="row")
            if i + 4 > len(self.navigation):
                page.div(class_='col-sm-3')
            else:
                page.div(class_='col-sm-3 mb-10')
            page.a(key, href=url, rel='external', target='_blank',
                   class_='btn btn-success btn-lg btn-block')
            page.div.close()
            if i % 4 == 3 or i == len(self.navigation) - 1:
                page.div.close()

        if self.plots:
            page.hr(class_='row-divider')
            page.div(class_='card border-light card-body scaffold shadow-sm')
            # parse frequencies
            freqs = [list(map(float, os.path.basename(p.href).split('_')[1:3]))
                     for p in self.plots[::2]]
            # build a grid of buttons as frequency markers
            ncols = 5
            linksize = 12 // ncols
            offset = 1
            i = 0
            page.p("Select a frequency range to jump to those plots:")
            for i, (f1, f2) in enumerate(freqs):
                if i % ncols == 0:
                    page.div(class_='row', style="margin-bottom: 10px;")
                if offset and i % ncols == 0:
                    page.div(class_="col-md-%d offset-md-%d"
                                    % (linksize, offset))
                else:
                    page.div(class_="col-md-%d" % linksize)
                page.div(class_="btn-group btn-group-justified")
                page.a('%s-%s Hz' % (f1, f2), href="javascript:;",
                       onclick="document.location.hash='%s-%s'" % (f1, f2),
                       role="button", type="button", class_="btn btn-warning")
                page.div.close()
                page.div.close()
                if i % ncols == ncols - 1:
                    page.div.close()
            if i % ncols != ncols - 1:
                page.div.close()
            page.hr(class_='row-divider')
            # reverse frequency order
            spectrograms = self.plots[::2]
            spectra = self.plots[1::2]
            for pair in list(zip(spectrograms, spectra, freqs))[::-1]:
                f = pair[-1]
                page.div(class_='row', id="%s-%s" % (f[0], f[1]))
                for p in pair[:2]:
                    page.div(class_="col-md-6")
                    page.a(href=p.href, class_="fancybox plot",
                           **{'data-fancybox-group': 1})
                    page.img(class_='img-fluid w-100', src=p.href)
                    page.a.close()
                    page.div.close()
                page.div.close()
            page.div.close()

        # write to file
        idx = self.states.index(state)
        with open(self.frames[idx], 'w') as fobj:
            fobj.write(str(page))
        return self.frames[idx]


register_tab(FscanTab)
