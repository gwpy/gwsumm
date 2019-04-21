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

from gwdetchar.utils import parse_html

from .. import html5

# global variables

URL = 'https://github.com/gwpy/gwsumm'
LOAD = """<script>
    $.ajax({
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


# test utilities

def test_load_state():
    state = html5.load_state('test')
    assert parse_html(str(state)) == parse_html(
        '<script>\nif (location.hash.length <= 1) {\n'
        '    $("#state_test").load_state("test");\n}\n'
        '</script>')


def test_load():
    # test local url
    content = html5.load('test')
    assert parse_html(content) == parse_html(
        '<script>$("#main").load("test");</script>')
    # test non-local url
    success = '$(\"#main\").html(data);'
    errormsg = ('alert(\"Cannot load content from %r, use '
                'browser console to inspect failure.\");' % URL)
    content = html5.load(URL)
    assert parse_html(content) == parse_html(LOAD % (URL, success, errormsg))
    # test with non-string error argument
    success = '$(\"#main\").html(data);'
    error = 'Failed to load content from %r' % URL
    errormsg = ('$(\"#main\").html(\"<div class=\'alert alert-warning\'>'
                '<p>%s</p></div>\");' % error)
    content = html5.load(URL, error=1)
    assert parse_html(content) == parse_html(LOAD % (URL, success, errormsg))


def test_load_custom():
    error = 'Error'
    success = 'document.write(\"Success\")'
    errormsg = ('$(\"#main\").html(\"<div class=\'alert alert-warning\'>'
                '<p>%s</p></div>\");' % error)
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
