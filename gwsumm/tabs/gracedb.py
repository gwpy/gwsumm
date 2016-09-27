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

from gwpy.time import from_gps

from .registry import (get_tab, register_tab)

from .. import html
from ..config import (GWSummConfigParser, NoOptionError)
from ..utils import (re_quote, vprint)

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__all__ = ['GraceDbTab']

LABELS = {
    'ADVOK': 'success',
    'ADVNO': 'danger',
    'H1OPS': 'info',
    'H1OK': 'success',
    'H1NO': 'danger',
    'L1OPS': 'info',
    'L1OK': 'success',
    'L1NO': 'danger',
    'DQV': 'danger',
    'INJ': 'warning',
    'EM_READY': 'success',
    'PE_READY': 'success',
}


class GraceDbTab(get_tab('default')):
    """Custom tab displaying a summary of GraceDb results.
    """
    type = 'gracedb'

    def __init__(self, name, url='https://gracedb.ligo.org',
                 query='External', columns=['gpstime', 'date', 'pipeline'],
                 headers=['GPS time', 'UTC time', 'Source'], rank='gpstime',
                 **kwargs):
        super(GraceDbTab, self).__init__(name, **kwargs)
        self.url = url
        self.query = query
        self.events = dict()
        self.headers = headers
        self.columns = columns
        self.rank = rank

    @classmethod
    def from_ini(cls, config, section, **kwargs):
        """Define a new `GraceDbTab` from a `ConfigParser`.
        """
        for key in ['url', 'query', 'rank']:
            try:
                kwargs.setdefault(
                    key, re_quote.sub('', config.get(section, key)))
            except NoOptionError:
                pass
        for key in ['columns', 'headers']:
            try:
                raw = config.get(section, key)
                l = eval(raw)
            except NoOptionError:
                continue
            except (SyntaxError, NameError, TypeError):
                l = [x.strip().rstrip() for x in raw.split(',')]
            kwargs.setdefault(key, l)
        return super(GraceDbTab, cls).from_ini(config, section, **kwargs)

    def process(self, config=GWSummConfigParser(), **kwargs):
        try:
            from ligo.gracedb.rest import GraceDb
        except ImportError as e:
            e.args = ('%s, this module is required to generate a GraceDbTab'
                      % str(e),)
            raise
        # query gracedb
        service_url = '%s/api/' % self.url
        connection = GraceDb(service_url=service_url)
        vprint("Connected to gracedb at %s\n" % connection.service_url)
        querystr = '%s %d .. %d' % (self.query, self.start, self.end)
        self.events[None] = list(connection.events(querystr))
        vprint("Recovered %d events for query %r\n"
               % (len(self.events[None]), querystr))
        if 'labels' in self.columns:
            for e in self.events[None]:
                e['labels'] = ', '.join(connection.event(
                    e['graceid']).json()['labels'])
            vprint("Downloaded labels\n")
        return super(GraceDbTab, self).process(config=config, **kwargs)

    def process_state(self, state, **kwargs):
        def in_state(event):
            return int(event['gpstime']) in state.active
        self.events[str(state)] = filter(in_state, self.events[None])
        reverse = self.rank not in ['gpstime', 'far']
        self.events[str(state)].sort(key=lambda x: x[self.rank],
                                     reverse=reverse)
        vprint("    Selected %d events\n" % len(self.events[str(state)]))

    def write_state_html(self, state):
        """Write the '#main' HTML content for this `GraceDbTab`.
        """
        page = html.markup.page()
        # build table of events
        page.div(class_='scaffold well')
        page.table(class_='table table-condensed table-hover table-striped')
        # thead
        page.thead()
        page.tr()
        for head in self.headers:
            page.th(head)
        page.tr.close()
        page.thead.close()
        # tbody
        page.tbody()
        for event in self.events[str(state)]:
            context = None
            try:
                l = event['labels'].split(', ')
            except (AttributeError, KeyError):
                pass
            else:
                for label in ['ADVNO', 'H1NO', 'L1NO', 'DQV', 'INJ',
                              'EM_READY']:
                    if label in l:
                        context = LABELS[label]
                        break
            if context is not None:
                page.tr(class_=context)
            for col in self.columns:
                if col == 'date':
                    page.td(from_gps(event['gpstime']).strftime(
                        '%B %d %Y, %H:%M:%S.%f')[:-3])
                    continue
                try:
                    v = event[col]
                except KeyError:
                    try:
                        v = event['extra_attributes']['GRB'][col]
                        assert v is not None
                    except (KeyError, AssertionError):
                        page.td('-')
                        continue
                if col == 'graceid':
                    page.td()
                    href='%s/events/view/%s' % (self.url, v)
                    page.a(v, href=href, target='_blank', rel='external')
                    page.td.close()
                elif col != 'gpstime' and isinstance(v, float):
                    page.td('%.3g' % v)
                else:
                    page.td(str(v))
            page.tr.close()
        page.tbody.close()
        page.table.close()
        if len(self.events[str(state)]) == 0:
            page.p("No events were recovered for this state.")
        page.div.close()  # scaffold well

        # query doc
        page.p("The above table was generated from a query to %s with the "
               "form <code>%s</code>." % (self.url, self.query))

        # reference the labelling
        page.h4("Labelling reference")
        page.p("Events in the above table may have a context based on "
               "its labels as follows:")
        contexts = set(LABELS.values())
        for c in contexts:
            labels = [k for k, v in LABELS.items() if v == c]
            labstr = ', '.join(['<b>%s</b>' % l for l in labels])
            page.p(labstr, class_='bg-%s' % c, style='width: auto;')

        # write to file
        idx = self.states.index(state)
        with open(self.frames[idx], 'w') as fobj:
            fobj.write(str(page))
        return self.frames[idx]

register_tab(GraceDbTab)
