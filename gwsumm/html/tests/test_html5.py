# -*- coding: utf-8 -*-
# Copyright (C) Alex Urban (2019)
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

"""Unit tests for gwsumm.html.html5
"""

__author__ = 'Alex Urban <alexander.urban@ligo.org>'

import os
import shutil

from markdown import markdown

from gwdetchar.utils import parse_html

from .. import html5

# global variables

URL = 'https://github.com/gwpy/gwsumm'
LOAD = """<script>
    jQuery.ajax({
        url : '%s',
        type : 'GET',
        success: function(data, statusText, jqhxr){%s},
        error: function(xhr, status, error){%s}
        });
</script>"""

BOX = """<div id="disqus_thread">
<script type="text/javascript">
    var disqus_shortname = "Test";
    var disqus_identifier = "test";
    var disqus_title = "Test";
    var disqus_url = "%s";

        (function() {
            var dsq = document.createElement('script');
            dsq.type = 'text/javascript'; dsq.async = true;
            dsq.src = '//' + disqus_shortname + '.disqus.com/embed.js';
            (document.getElementsByTagName('head')[0] ||
             document.getElementsByTagName('body')[0]).appendChild(dsq);
        })();

</script>
<noscript>Please enable JavaScript to view the</noscript>
<a href="https://disqus.com/?ref_noscript">comments powered by Disqus</a>
</div>"""

CONTENTS = """Heading
=======

This is a test.

* Bullet 1
* Bullet 2"""

DIALOG = ('<button title="Test" id="id-btn" class="btn-float btn-open shadow" '
          'data-id="#id">T</button>\n<div class="dialog" title="Test" '
          'id="id">\n%s\n</div>') % markdown(CONTENTS)

OVERLAY = (
    '<button title="Overlay figures" id="overlay-btn" class="btn-float '
    'btn-open shadow" data-id="#overlay"><i class="fas fa-layer-group"></i>'
    '</button>\n<div class="dialog" title="Overlay figures" id="overlay">\n'
    '<h1>Overlay figures for easy comparison</h1>\n<hr class="row-divider" />'
    '\n<div class="row" id="overlay-outer">\n<div class="col-md-4">\n'
    '<div class="card card-body shadow-sm" id="overlay-info">\n'
    '<h4>Instructions</h4>\n%s\n</div>\n</div>\n<div class="col-md-8">\n'
    '<div class="text-center">\n<a title="Overlay all selected figures" '
    'class="btn btn-light shadow-sm" id="overlay-figures">Overlay</a>\n'
    '<a title="Download overlay figure" class="btn btn-light shadow-sm" '
    'id="download-overlay">Download</a>\n<a title="Clear all figure '
    'selections" class="btn btn-light shadow-sm" id="clear-figures">Clear</a>'
    '\n</div>\n<br />\n<canvas id="overlay-canvas" />\n</div>\n</div>\n'
    '</div>') % markdown(html5.OVERLAY_INSTRUCTIONS)


# test utilities

def test_expand_path():
    assert html5._expand_path(URL) == URL


def test_load_state():
    state = html5.load_state('test')
    assert parse_html(str(state)) == parse_html(
        '<script>\nif (location.hash.length <= 1) {\n'
        '    jQuery("#state_test").load_state("test");\n}\n'
        '</script>')


def test_load():
    # test local url
    content = html5.load('test')
    assert parse_html(content) == parse_html(
        '<script>jQuery("#main").load("test");</script>')
    # test non-local url
    success = 'jQuery(\"#main\").html(data);'
    errormsg = ('alert(\"Cannot load content from %r, use '
                'browser console to inspect failure.\");' % URL)
    content = html5.load(URL)
    assert parse_html(content) == parse_html(LOAD % (URL, success, errormsg))
    # test with non-string error argument
    success = 'jQuery(\"#main\").html(data);'
    error = 'Failed to load content from %r' % URL
    errormsg = ('jQuery(\"#main\").html(\"<div class=\'alert alert-warning '
                'text-justify shadow-sm\'><p>%s</p></div>\");' % error)
    content = html5.load(URL, error=1)
    assert parse_html(content) == parse_html(LOAD % (URL, success, errormsg))


def test_load_custom():
    error = 'Error'
    success = 'document.write(\"Success\")'
    errormsg = ('jQuery(\"#main\").html(\"<div class=\'alert alert-warning '
                'text-justify shadow-sm\'><p>%s</p></div>\");' % error)
    # test local url
    content = html5.load('test', success=success, error=error)
    assert parse_html(content) == parse_html(
        LOAD % ('test', success, errormsg))
    # test non-local url
    content = html5.load(URL, success=success, error=error)
    assert parse_html(content) == parse_html(LOAD % (URL, success, errormsg))


def test_comments_box():
    box = html5.comments_box('Test', identifier='test', title='Test', url=URL)
    assert parse_html(str(box)) == parse_html(BOX % URL)


def test_ldvw_qscan_single():
    button = html5.ldvw_qscan('X1:TEST', 0)
    assert parse_html(str(button)) == parse_html(
        '<a href="https://ldvw.ligo.caltech.edu/ldvw/view?act=doplot&amp;'
        'chanName=X1:TEST&amp;doQxfrm=yes&amp;strtTime=0&amp;&amp;'
        'qxfrm_pltfrq=10 inf&amp;qxfrm_srchqrng=4 100&amp;'
        'qxfrm_plttimes=0.5 2 8" target="_blank" rel="external" '
        'class="btn btn-outline-secondary btn-sm" title="Q-scan X1:TEST at 0 '
        'via LDVW">Q-scan</a>')


def test_ldvw_qscan_batch():
    button = html5.ldvw_qscan('X1:TEST', (0,))
    assert parse_html(str(button)) == parse_html(
        '<a href="https://ldvw.ligo.caltech.edu/ldvw/Wdq?submitAct=go&amp;'
        'wdq_ifo=X1&amp;wdq_cmap=viridis&amp;wdq_gps=0&amp;wdq_prog=py-Omega&'
        'amp;goBtn=goBtn" target="_blank" rel="external" class="btn '
        'btn-outline-secondary btn-sm" title="Batch-process omega scans of '
        'the loudest triggers via LDVW">Launch omega scans</a>')


def test_dialog_box(tmpdir):
    mdfile = os.path.join(str(tmpdir), 'test.md')
    with open(mdfile, 'w') as f:
        f.write(CONTENTS)
    box = html5.dialog_box(mdfile, 'Test', 'id', 'T')
    assert parse_html(str(box)) == parse_html(DIALOG)
    shutil.rmtree(str(tmpdir), ignore_errors=True)


def test_overlay_canvas():
    box = html5.overlay_canvas()
    assert parse_html(str(box)) == parse_html(OVERLAY)
