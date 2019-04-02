# -*- coding: utf-8 -*-
# Copyright (C) Duncan Macleod (2016)
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

"""HTML <table> helpers
"""

from MarkupPy.markup import page

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'


def table(headers, data, caption=None, separator='', id=None, **class_):
    """Write a <table> with one row of headers and many rows of data

    Parameters
    ----------
    headers : `list`
        list of column header names

    data : `list` of `lists`
        list of column data lists, for ``m`` rows and ``n`` columns, this
        should have dimensions ``m x n``

    caption : `str`, optional
        content for this table's `<caption>`

    **class_
        class attribute declarations for each tag used in the table,
        any of `table`, `thead`, `tbody`, `tr`, `th`, `td`, `caption`

    Returns
    -------
    table : `markup.page`
        a formatted HTML page object containing the `<table>`
    """
    class_.setdefault('table',
                      'table table-hover table-condensed table-responsive')
    # unwrap class declarations (so we don't get empty class attributes)
    kwargs = {}
    for tag in ['table', 'thead', 'tbody', 'tr', 'th', 'td', 'caption']:
        try:
            kwargs[tag] = {'class_': class_.pop(tag)}
        except KeyError:
            kwargs[tag] = {}

    # create table and add caption
    p = page(separator=separator)
    if id is not None:
        kwargs['table']['id_'] = id
    p.table(**kwargs['table'])
    if caption:
        p.caption(caption, **kwargs['caption'])

    # write headers
    p.thead(**kwargs['thead'])
    p.tr(**kwargs['tr'])
    for th in headers:
        p.th(th, scope='col', **kwargs['th'])
    p.tr.close()
    p.thead.close()

    # write body
    p.tbody(**kwargs['tbody'])
    for row in data:
        p.tr(**kwargs['tr'])
        for td in row:
            p.td(td, **kwargs['td'])
        p.tr.close()
    p.tbody.close()

    p.table.close()

    if id is not None:
        p.button(
            'Export to CSV', class_='btn btn-default btn-table',
            onclick="exportTableToCSV('{name}.csv', '{name}')".format(
                name=id))
    return p
