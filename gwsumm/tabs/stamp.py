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
import re
import glob

from six import string_types

from .registry import (get_tab, register_tab)

from .. import html
from ..mode import Mode
from ..plot import get_plot
from ..config import GWSummConfigParser

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__all__ = ['StampPEMTab']

base = get_tab('default')
SummaryPlot = get_plot(None)


class StampPEMTab(base):
    """Custom tab displaying a summary of StampPEM results.
    """
    type = 'stamp'

    def __init__(self, *args, **kwargs):
        if kwargs['mode'] != Mode.day:
            raise RuntimeError("StampPEMTab is only available in %s mode."
                               % Mode.day.name)
        super(StampPEMTab, self).__init__(*args, **kwargs)

    @classmethod
    def from_ini(cls, config, section, **kwargs):
        """Define a new `StampPEMTab` from a `ConfigParser`.
        """
        # parse generic configuration
        new = super(StampPEMTab, cls).from_ini(config, section, **kwargs)
        new.set_layout([2])

        # work out day directory and url
        new.directory = os.path.normpath(config.get(section, 'base-directory'))
        try:
            home_, postbase = new.directory.split('/public_html/', 1)
        except ValueError as e:
            e.args = ('Stamp PEM directory not under \'public_html\', '
                      'cannot format linkable URL',)
            raise
        else:
            user = os.path.split(home_)[1]
            new.url = '/~%s/%s' % (user, postbase.rstrip('/'))
        return new

    def process(self, config=GWSummConfigParser(), **kwargs):
        # find all plots
        self.plots = []
        if isinstance(self.directory, string_types):
            plots = sorted(
                glob.glob(os.path.join(self.directory, 'DAY_*.png')),
                key=lambda p: float(re.split(r'[-_]', os.path.basename(p))[1]))
            for p in plots:
                pname = os.path.split(p)[1]
                self.plots.append(SummaryPlot(
                    src=os.path.join(self.url, pname),
                    href=os.path.join(self.url,
                                      pname.replace('.png', '.html'))))

    def write_state_html(self, state):
        """Write the '#main' HTML content for this `StampPEMTab`.
        """
        page = markup.page()

        a = markup.oneliner.a('analysis', href=self.url+'/',
                              class_='alert-link', rel='external',
                              target='_blank')
        if not(os.path.isdir(self.directory)):
            page.div(class_='alert alert-warning', role='alert')
            page.p("No %s was performed for this period, "
                   "please try again later." % a)
            page.p("If you believe these data should have been found, please "
                   "contact %s."
                   % markup.oneliner.a('the DetChar group',
                                       class_='alert-link',
                                       href='mailto:detchar@ligo.org'))
            page.div.close()

        elif not self.plots:
            page.div(class_='alert alert-warning', role='alert')
            page.p("This %s produced no plots." % a)
            page.p("If you believe these data should have been found, please "
                   "contact %s."
                   % markup.oneliner.a('the DetChar group',
                                       class_='alert-link',
                                       href='mailto:detchar@ligo.org'))
            page.div.close()

        else:
            page.add(str(self.scaffold_plots(
                aclass='fancybox-stamp plot',
                **{'data-fancybox-type': 'iframe'})))
            page.hr(class_='row-divider')

            # link full results
            page.hr(class_='row-divider')
            page.div(class_='btn-group')
            page.a('Click here for the full Stamp PEM results',
                   href=self.url+'/', rel='external', target='_blank',
                   class_='btn btn-default btn-info btn-xl')
            page.div.close()

        # write to file
        idx = self.states.index(state)
        with open(self.frames[idx], 'w') as fobj:
            fobj.write(str(page))
        return self.frames[idx]


register_tab(StampPEMTab)
