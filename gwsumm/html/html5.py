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

"""HTML5-specific extensions
"""

import re
import os.path

from pathlib import Path
from urllib.parse import urlparse

from MarkupPy import markup
from markdown import markdown

# global variables

DOCTYPE = '<!DOCTYPE html>'

COMMENTS_JS = """
        (function() {
            var dsq = document.createElement('script');
            dsq.type = 'text/javascript'; dsq.async = true;
            dsq.src = '//' + disqus_shortname + '.disqus.com/embed.js';
            (document.getElementsByTagName('head')[0] ||
             document.getElementsByTagName('body')[0]).appendChild(dsq);
        })();
"""

OVERLAY_INSTRUCTIONS = """\
This tool is designed to help visually correlate multiple types of data by
laying transparent figures atop one another. To use it, follow these steps:

* From the main page(s), drag-and-drop figures onto the button with this icon:
  <i class="fas fa-layer-group"></i>
* Expand this dialog box and use the `Overlay` button (to the right) to render
  a new figure by overlaying the ones selected
* Add more figures using the steps above, click `Overlay` to re-render, and
  use the `Download` button to save to a local file
* Use the `Clear` button to clear all figure selections

**Note:** figures are remembered across tabs, so multiple times and subsystems
can be compared at once.\
"""


# -- utilities ----------------------------------------------------------------

def _expand_path(path):
    """Expand a server path that may contain symbolic links
    """
    subbed = Path(re.sub(r'^/\~(.*?)/', r'/home/\1/public_html/', path))
    resolved = subbed.resolve() if subbed.exists() else subbed
    return re.sub(r'^/home/(.*?)/public_html/', r'/~\1/',
                  str(resolved) if resolved.exists() else path)


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
    page.add('    jQuery("#state_%s").load_state("%s");' % (id_, url))
    page.add('}')
    page.script.close()
    return page


def load(url, id_='main', error=False, success=None):
    """Construct the HTML script required to load a url into the
    HTML element with the given unique ``id_``.
    """
    ps = urlparse(url)
    if not ps.netloc and not error:
        return markup.given_oneliner.script('jQuery("#%s").load("%s");'
                                            % (id_, url))
    elif ps.netloc and not error:
        error = ('alert("Cannot load content from %r, use browser console '
                 'to inspect failure.");' % url)
    else:
        if not isinstance(error, (str, markup.page)):
            error = 'Failed to load content from %r' % url
        error = ('jQuery("#%s").html("<div class=\'alert alert-warning '
                 'text-justify shadow-sm\'><p>%s</p></div>");' % (id_, error))
    if success is None:
        success = 'jQuery("#%s").html(data);' % id_
    url = _expand_path(url)
    return markup.given_oneliner.script("""
    jQuery.ajax({
        url : '%s',
        type : 'GET',
        success: function(data, statusText, jqhxr){%s},
        error: function(xhr, status, error){%s}
        });\n""" % (url, success, error))


def comments_box(name, identifier=None, title=None, url=None):
    """Generate a Disqus comments box
    """
    page = markup.page()
    page.div(id_='disqus_thread')
    page.script(type='text/javascript')
    page.add('    var disqus_shortname = "%s";' % name)
    if identifier:
        page.add('    var disqus_identifier = "%s";' % identifier)
    if title:
        page.add('    var disqus_title = "%s";' % title)
    if url:
        page.add('    var disqus_url = "%s";' % url)
    page.add(COMMENTS_JS)
    page.script.close()
    page.noscript("Please enable JavaScript to view the")
    page.a("comments powered by Disqus",
           href="https://disqus.com/?ref_noscript")
    page.div.close()
    return page


def ldvw_qscan(channel, time, fmin=10, fmax='inf', qmin=4, qmax=100):
    """Generate a Q-scan through LIGO DataViewer Web (LDVW)
    """
    channel = str(channel)
    if isinstance(time, (tuple, list)):
        label = 'Launch omega scans'
        title = 'Batch-process omega scans of the loudest triggers via LDVW'
        times = '&'.join('wdq_gps=' + str(t) for t in time)
        query = ('Wdq?submitAct=go&wdq_ifo=%s&wdq_cmap=viridis&%s&'
                 'wdq_prog=py-Omega&goBtn=goBtn') % (channel[:2], times)
    else:
        label = 'Q-scan'
        title = 'Q-scan {0} at {1} via LDVW'.format(channel, time)
        query = ('view?act=doplot&chanName={0}&doQxfrm=yes&strtTime={1}&'
                 '&qxfrm_pltfrq={2} {3}&qxfrm_srchqrng={4} {5}&'
                 'qxfrm_plttimes=0.5 2 8').format(
                     channel, time, fmin, fmax, qmin, qmax)
    uri = 'https://ldvw.ligo.caltech.edu/ldvw/{0}'.format(query)
    return markup.oneliner.a(
        label, href=uri, target='_blank', rel='external',
        class_='btn btn-outline-secondary btn-sm', title=title)


def dialog_box(content, title, id_, btntxt):
    """Generate a dialog box to be loaded modal atop the main page

    Parameters
    ----------
    content : `str`
        either raw markdown text or the path to a file containing markdown,
        this will be rendered in HTML as the contents of the dialog box

    title : `str`
        title to display atop the dialog box

    id_ : `str`
        unique identifier for the dialog box

    btntxt : `str`
        text (usually a single character) to appear inside a sticky button
        that opens the dialog box

    Returns
    -------
    page : `~MarkupPy.markup.page`
        fully rendered HTML containing the dialog box
    """
    btnargs = {
        'title': title,
        'id_': '-'.join([id_, 'btn']),
        'class_': 'btn-float btn-open shadow',
        'data-id': '#' + id_,
    }
    page = markup.page()
    page.button(btntxt, **btnargs)
    page.div(
        class_='dialog',
        title=title,
        id_=id_,
    )
    if os.path.isfile(content):
        with open(content, 'r') as source:
            content = source.read()
    page.add(markdown(str(content)))
    page.div.close()
    return page


def overlay_canvas():
    """Generate a dialog box allowing users to select and overlay plots

    Returns
    -------
    page : `~MarkupPy.markup.page`
        fully rendered HTML containing the dialog box
    """
    page = markup.page()
    page.h1('Overlay figures for easy comparison')
    page.hr(class_='row-divider')
    page.div(class_='row', id_='overlay-outer')
    page.div(class_='col-md-4')
    page.div(class_='card card-body shadow-sm', id_='overlay-info')
    page.h4('Instructions')
    page.add(markdown(OVERLAY_INSTRUCTIONS))
    page.div.close()  # card card-body shadow-sm
    page.div.close()  # col-md-4
    page.div(class_='col-md-8')
    page.div(class_='text-center')
    page.a('Overlay', title='Overlay all selected figures',
           class_='btn btn-light shadow-sm', id_='overlay-figures')
    page.a('Download', title='Download overlay figure',
           class_='btn btn-light shadow-sm', id_='download-overlay')
    page.a('Clear', title='Clear all figure selections',
           class_='btn btn-light shadow-sm', id_='clear-figures')
    page.div.close()  # text-center
    page.br()
    page.add(markup.oneliner.canvas(id_='overlay-canvas'))
    page.div.close()  # col-md-8
    page.div.close()  # row
    return dialog_box(
        str(page), title='Overlay figures', id_='overlay',
        btntxt=markup.oneliner.i('', class_='fas fa-layer-group'))
