# -*- coding: utf-8
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

"""HTML5 specific extensions
"""

import os.path
from urlparse import urlparse

from . import markup
from ..utils import re_cchar

DOCTYPE = '<!DOCTYPE html>'


def load_state(url):
    """Construct the HTML script required to load the Tab HTML
    for a given :class:`~gwsumm.state.core.SummaryState`

    Parameters
    ----------
    url : `str`
        path (relative to <base>) of HTML to load
    id_ : `str`, optional, default: '#main'
        <div> 'id' in which to load HTML

    Returns
    -------
    HTML : `str`
        HTML one-liner with script loading
    """
    id_ = os.path.splitext(os.path.basename(url))[0]
    page = markup.page()
    page.script()
    page.add('if (location.hash.length <= 1) {')
    page.add('    $("#state_%s").load_state("%s");' % (id_, url))
    page.add('}')
    page.script.close()
    return page


def load(url, id_='main'):
    """Construct the HTML script required to load a url into the
    HTML element with the given unique ``id_``.
    """
    ps = urlparse(url)
    if ps.netloc:
        return markup.given_oneliner.script("""
    $.ajax({
        url : %r,
        type : 'GET',
        success: function(data){$("#%s").html(data);},
        error: function(xhr, status, error){
                   alert("Cannot load content from %r, use browser console" +
                         " to inspect failure.");
                   }
        });\n""" % (url, id_, url))
    else:
        return markup.given_oneliner.script('$("#%s").load("%s");' % (id_, url))
