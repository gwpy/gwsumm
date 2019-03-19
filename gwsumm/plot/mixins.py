# -*- coding: utf-8 -*-
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

"""
"""

import abc
import os.path
import re

from six.moves import StringIO

from lxml import etree

from .utils import usetex_tex

re_bit_label = re.compile(r'\[(?P<idx>.*)\] (?P<label>.*)')
re_source_label = re.compile(r'(?P<label>.*) \[(?P<flag>.*)\]')

HOVERSCRIPT = """
<script type="text/ecmascript">
<![CDATA[    function init(evt) {
      if ( window.svgDocument == null ) {
        svgDocument = evt.target.ownerDocument;
      }
      var labels = svgDocument.getElementsByClassName('mpl-label');
      for (var i=0; labels[i]; i++) {
        labels[i].setAttribute('visibility', 'hidden');
      }
    }

    function ShowLabel(obj) {
      var cur = obj.id.substr(obj.id.lastIndexOf('_')+1);
      var tip = svgDocument.getElementById('label_' + cur);
      tip.setAttribute('visibility',"visible")
    }

    function HideLabel(obj) {
      var cur = obj.id.substr(obj.id.lastIndexOf('_')+1);
      var tip = svgDocument.getElementById('label_' + cur);
      tip.setAttribute('visibility',"hidden")
    }]]>
  </script>

"""

HTML_WRAPPER = """
<html>
<body>
  <style>
    body { margin: 0 !important }
  </style>
</head>
<div>
<object type="image/svg+xml" data="%s"
Your browser cannot display this SVG
</object>
</div>
</body>
</html>
"""


class SvgMixin(object):
    __metaclass__ = abc.ABCMeta

    def __init__(self, *args, **kwargs):
        super(SvgMixin, self).__init__(*args, **kwargs)
        self.preview_labels = False

    def finalize(self, outputfile=None, close=True, **savekwargs):
        if outputfile is None:
            outputfile = self.outputfile
        # make SVG
        if outputfile.endswith('.svg'):
            self.draw_svg(outputfile)
        # or continue as normal
        else:
            return super(SvgMixin, self).finalize(
                outputfile=outputfile, close=close, **savekwargs)

    @abc.abstractmethod
    def draw_svg(self, outputfile):
        pass

    def finalize_svg(self, tree, outputfile, script=None):
        if script:
            tree.insert(0, etree.XML(script))
        etree.ElementTree(tree).write(outputfile)
        # write HTML wrapper
        html = self.outputfile.replace('.svg', '.html')
        with open(html, 'w') as f:
            f.write(HTML_WRAPPER % os.path.basename(self.outputfile))
        return html


class DataLabelSvgMixin(SvgMixin):
    def draw_svg(self, outputfile):
        for ax in self.plot.axes:
            for line in ax.lines:
                line.set_rasterized(True)

        # render image
        super(SvgMixin, self).finalize(
            outputfile=outputfile.replace('.svg', '.png'), close=False)

        # make new text labels for the channel names
        ax = self.plot.axes[0]
        leg = ax.legend_
        texts = []
        if leg is not None:
            for i, (text, line) in enumerate(
                    zip(leg.get_texts(), leg.get_lines())):
                try:
                    label, source = re_source_label.match(
                        text.get_text()).groups()
                except (AttributeError, ValueError):
                    continue
                channel = usetex_tex(str(source))
                text.set_text(label)
                t2 = ax.text(
                    0.994, 1.02, channel, ha='right', va='bottom',
                    fontsize=text.get_fontsize(), zorder=1000,
                    transform=ax.transAxes,
                    bbox={'facecolor': 'white', 'edgecolor': 'lightgray',
                          'pad': 10.})
                text.set_gid('leg_text_%d' % i)
                line.set_gid('leg_patch_%d' % i)
                t2.set_gid('label_%d' % i)
                texts.append(t2)

        # tmp save
        f = StringIO()
        self.plot.save(f, format='svg')

        # parse svg
        tree, xmlid = etree.XMLID(f.getvalue())
        tree.set('onload', 'init(evt)')

        # add effects
        for i in range(len(texts)):
            pl = xmlid['leg_patch_%d' % i]
            ll = xmlid['leg_text_%d' % i]
            tl = xmlid['label_%d' % i]
            pl.set('cursor', 'pointer')
            pl.set('onmouseover', "ShowLabel(this)")
            pl.set('onmouseout', "HideLabel(this)")
            ll.set('cursor', 'pointer')
            ll.set('onmouseover', "ShowLabel(this)")
            ll.set('onmouseout', "HideLabel(this)")
            tl.set('class', 'mpl-label')
            tl.set('visibility', 'hidden')

        return self.finalize_svg(tree, outputfile, script=HOVERSCRIPT)


class SegmentLabelSvgMixin(SvgMixin):
    def draw_svg(self, outputfile):
        # render image
        super(SvgMixin, self).finalize(outputfile=StringIO(), close=False)
        ax = self.plot.axes[0]
        collections = [c for c in ax.collections if hasattr(c, '_ignore')]
        # reset labels
        labels = {}
        for i, t in enumerate(ax.get_yaxis().get_ticklabels()):
            text = t.get_text()
            if not text:
                continue
            x, y = t.get_position()
            m1 = re_bit_label.match(text)
            m2 = re_source_label.match(text)
            if m1:
                idx, label = m1.groups()
                t2 = ax.text(0.01, y, label, ha='left',
                             fontsize=t.get_fontsize(), va='center',
                             transform=t.get_transform())
                t2.set_bbox({'alpha': 0.5, 'facecolor': 'white',
                             'edgecolor': 'none'})
                t.set_text(idx)
                t.set_bbox(None)
                t.set_position((0, y))
                t.set_ha('right')
                t.set_fontsize('14')
                labels[text] = (idx, t2)
            elif m2:
                label, flag = m2.groups()
                t2 = ax.text(x, y, text, ha='left', bbox=t._bbox.copy(),
                             fontsize=t.get_fontsize(), va='center',
                             transform=t.get_transform())
                t.set_text(label)
                labels[text] = (label, t2)
        ax._insetlabels = None

        j = 0
        for i, collection in enumerate(collections):
            try:
                idx, tickl = labels[collection.get_label()]
            except KeyError:
                continue
            else:
                labels[collection] = tickl
            collection.set_label(idx)
            if not tickl.get_gid():
                tickl.set_gid('label_%d' % j)
                j += 1
            collection.set_gid('collection_%d_%s' % (i, tickl.get_gid()))

        # render as and parse SVG
        f = StringIO()
        self.plot.save(f, format='svg')
        tree, xmlid = etree.XMLID(f.getvalue())
        tree.set('onload', 'init(evt)')
        # set mouse events for visible labels
        for i, collection in enumerate(collections):
            try:
                tickl = labels[collection]
            except KeyError:
                continue
            cel = xmlid[collection.get_gid()]
            cel.set('cursor', 'pointer')
            cel.set('onmouseover', "ShowLabel(this)")
            cel.set('onmouseout', "HideLabel(this)")
            tel = xmlid[tickl.get_gid()]
            tel.set('class', 'mpl-label')
            if not self.preview_labels:
                tel.set('visibility', 'hidden')

        # ignore hover on text events
        for key in xmlid:
            if key.startswith('text_') or key.startswith('label_'):
                xmlid[key].set('style', 'pointer-events: none;')

        return self.finalize_svg(tree, outputfile, script=HOVERSCRIPT)
