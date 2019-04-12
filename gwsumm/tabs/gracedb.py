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

from collections import OrderedDict
from configparser import NoOptionError

from MarkupPy import markup

from gwpy.time import from_gps

from .registry import (get_tab, register_tab)
from .. import html
from ..config import GWSummConfigParser
from ..utils import (re_quote, vprint)

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__all__ = ['GraceDbTab']

LABELS = OrderedDict()
LABELS["danger"] = {
    "ADVNO",
    "DQV",
    "H1NO",
    "L1NO",
    "V1NO",
}
LABELS["warning"] = {
    "INJ",
}
LABELS["success"] = {
    'GCN_PRELIM_SENT',
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
        self.query = "{} {} .. {}".format(
            query,
            int(self.start),
            int(self.end),
        )
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
                val = eval(raw)
            except NoOptionError:
                continue
            except (SyntaxError, NameError, TypeError):
                val = [x.strip().rstrip() for x in raw.split(',')]
            kwargs.setdefault(key, val)
        return super(GraceDbTab, cls).from_ini(config, section, **kwargs)

    def process(self, config=GWSummConfigParser(), **kwargs):
        try:
            from ligo.gracedb.rest import GraceDb
            from ligo.gracedb.exceptions import HTTPError
        except ImportError as e:
            e.args = ('%s, this module is required to generate a GraceDbTab'
                      % str(e),)
            raise
        # query gracedb
        service_url = '%s/api/' % self.url
        connection = GraceDb(service_url=service_url)
        vprint("Connected to gracedb at %s\n" % service_url)
        try:
            self.events[None] = list(connection.superevents(self.query))
            self._query_type = "S"
        except HTTPError:
            self.events[None] = list(connection.events(self.query))
            event_method = connection.event
            eventid_name = "graceid"
            self._query_type = "E"
        else:
            event_method = connection.superevent
            eventid_name = "superevent_id"
            for event in self.events[None]:  # get preferred event parameters
                event.update(connection.event(
                    event["preferred_event"],
                ).json())
        vprint("Recovered %d events for query %r\n"
               % (len(self.events[None]), self.query))
        if 'labels' in self.columns:
            for e in self.events[None]:
                e['labels'] = ', '.join(event_method(
                    e[eventid_name]).json()['labels'])
            vprint("Downloaded labels\n")
        return super(GraceDbTab, self).process(config=config, **kwargs)

    def process_state(self, state, **kwargs):
        def in_state(event):
            return int(event['gpstime']) in state.active
        self.events[str(state)] = list(filter(in_state, self.events[None]))
        reverse = self.rank not in ['gpstime', 'far']
        self.events[str(state)].sort(key=lambda x: x[self.rank],
                                     reverse=reverse)
        vprint("    Selected %d events\n" % len(self.events[str(state)]))

    def write_state_html(self, state):
        """Write the '#main' HTML content for this `GraceDbTab`.
        """
        page = markup.page()
        # build table of events
        page.div(class_='scaffold well')
        page.table(class_='table table-condensed table-hover table-striped',
                   id_='gracedb')
        # thead
        page.thead()
        page.tr()
        for head in self.headers:
            page.th(head)
        page.tr.close()
        page.thead.close()
        # tbody
        page.tbody()
        for event in sorted(self.events[str(state)],
                            key=lambda e: e['gpstime']):
            context = None
            try:
                labs = set(event['labels'].split(', '))
            except (AttributeError, KeyError):
                pass
            else:
                for ctx, labels in LABELS.items():
                    if (
                            ctx == "success" and labs.union(labels) == labs or
                            labs.intersection(labels)
                    ):
                        context = ctx
                        break
            if context:
                page.tr(class_=context)
            else:
                page.tr()
            for col in self.columns:
                if col == 'date':
                    gpskey = 't_0' if 'superevent_id' in event else 'gpstime'
                    page.td(from_gps(event[gpskey]).strftime(
                        '%B %d %Y %H:%M:%S.%f',
                    )[:-3])
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
                if col in ("graceid", "superevent_id", "preferred_event"):
                    page.td()
                    tag = "superevents" if col == "superevent_id" else "events"
                    href = '{}/{}/view/{}'.format(self.url, tag, v)
                    title = "GraceDB {} page for {}".format(tag[:-1], v)
                    page.a(v, title=title, href=href, target='_blank',
                           rel='external', class_="btn btn-info btn-xs")
                    page.td.close()
                elif col not in ("gpstime", "t_0") and isinstance(v, float):
                    page.td('%.3g' % v)
                elif col == "labels":
                    page.td(", ".join(['<samp>%s</samp>' % l for l in sorted(labs)]))
                else:
                    page.td(str(v))
            page.tr.close()
        page.tbody.close()
        page.table.close()
        if len(self.events[str(state)]) == 0:
            page.p("No events were recovered for this state.")
        else:
            page.button(
                'Export to CSV', class_='btn btn-default btn-table',
                onclick="exportTableToCSV('{name}.csv', '{name}')".format(
                    name='gracedb'))
        page.div.close()  # scaffold well

        # query doc
        qurl = "{}/search/?query={}&query_type={}&results_format=S".format(
            self.url,
            self.query.replace(" ", "+"),
            getattr(self, "_query_type", "E"),
        )
        qlink = markup.oneliner.a(
            "here",
            href=qurl,
            target="_blank",
        )
        page.p("The above table was generated from a query to {} with the "
               "form <code>{}</code>. To view the results of the same query "
               "via the GraceDB web interface, click {}.".format(
                   self.url, self.query, qlink),
        )

        # reference the labelling
        page.h4("Labelling reference")
        page.p("Events in the above table may have a context based on "
               "its labels as follows:")
        for c, labels in LABELS.items():
            labstr = ', '.join(['<samp>%s</samp>' % l for l in sorted(labels)])
            page.p(labstr, class_='bg-%s' % c, style='width: auto;')

        # write to file
        idx = self.states.index(state)
        with open(self.frames[idx], 'w') as fobj:
            fobj.write(str(page))
        return self.frames[idx]


register_tab(GraceDbTab)
