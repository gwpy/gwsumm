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

from urllib.parse import urlparse

from MarkupPy import markup

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


def load(url, id_='main', error=False, success=None):
    """Construct the HTML script required to load a url into the
    HTML element with the given unique ``id_``.
    """
    ps = urlparse(url)
    if not ps.netloc and not error:
        return markup.given_oneliner.script('$("#%s").load("%s");'
                                            % (id_, url))
    elif ps.netloc and not error:
        error = ('alert("Cannot load content from %r, use browser console '
                 'to inspect failure.");' % url)
    else:
        if not isinstance(error, (str, markup.page)):
            error = 'Failed to load content from %r' % url
        error = ('$("#%s").html("<div class=\'alert alert-warning\'>'
                 '<p>%s</p></div>");' % (id_, error))
    if success is None:
        success = '$("#%s").html(data);' % id_
    return markup.given_oneliner.script("""
    $.ajax({
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
    page.add("""
        (function() {
            var dsq = document.createElement('script');
            dsq.type = 'text/javascript'; dsq.async = true;
            dsq.src = '//' + disqus_shortname + '.disqus.com/embed.js';
            (document.getElementsByTagName('head')[0] ||
             document.getElementsByTagName('body')[0]).appendChild(dsq);
        })();
    """)
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
        class_='btn btn-default btn-xs', title=title)


def dialog_box(content, title, id_):
    """Generate a dialog box to be loaded modal atop the main page
    """
    page = markup.page()
    page.add('<dialog id="%s">' % id_)  # MarkupPy does not support dialog
    page.a('&#x2715;', title='Close', onclick="closeDialog('%s')" % id_,
           class_='btn btn-default pull-right', **{'aria-label': 'Close'})
    page.h1(title)
    page.add('<hr class="row-divider">')
    if isinstance(content, str):
        content = markup.oneliner.p(content) if (
            not content.startswith('<')) else content
    page.add(str(content))
    page.add('</dialog>')
    return page
