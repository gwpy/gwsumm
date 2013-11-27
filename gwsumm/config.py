
"""Thin wrapper around configparser
"""

try:
    from configparser import *
except ImportError:
    from ConfigParser import *
    from ConfigParser import __all__
else:
    from configparser import __all__

__all__.append('GWSummConfigParser')


class GWSummConfigParser(ConfigParser):
    def nditems(self, section):
        items = ConfigParser.items(self, section)
        return [i for i in items if i[0] not in self._defaults]
