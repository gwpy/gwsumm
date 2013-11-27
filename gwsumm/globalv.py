
"""Set of global memory variables for GWSumm package
"""

import time

from gwpy.segments import DataQualityDict

CHANNELS = {}
STATES = {}

DATA = {}
SPECTROGRAMS = {}
SPECTRUM = {}
SEGMENTS = DataQualityDict()

VERBOSE = False
PROFILE = False
START = time.time()

# run time variables
MODE = None
WRITTEN_PLOTS = []