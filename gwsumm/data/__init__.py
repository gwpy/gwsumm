# -*- coding: utf-8 -*-
# Copyright (C) Duncan Macleod (2013-2016)
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

"""Methods and classes for loading and pre-processing data

Each of the sub-modules are designed to read or create the data requested
only once, with the containers from the `globalv` module used as a storage
buffer for each unique data type
"""

# read TimeSeries data
from .timeseries import (
    _urlpath,
    _get_timeseries_dict,
    sieve_cache,
    find_frames,
    find_best_frames,
    find_frame_type,
    frame_trend_type,
    get_channel_type,
    exclude_short_trend_segments,
    all_adc,
    get_timeseries_dict,
    locate_data,
    get_timeseries,
    add_timeseries,
    resample_timeseries_dict,
    filter_timeseries,
)

# generate Spectrograms and FrequencySeries
from .spectral import (
    _get_spectrum,
    _get_spectrogram,
    get_spectrogram,
    add_spectrogram,
    get_spectrograms,
    size_for_spectrogram,
    apply_transfer_function_series,
    get_spectrum,
)

# generate coherence Spectrograms and Spectra
from .coherence import (
    _get_from_list,
    _get_coherence_spectrogram,
    get_coherence_spectrogram,
    get_coherence_spectrum,
    add_coherence_component_spectrogram,
    get_coherence_spectrograms,
    complex_percentile,
)

# generate TimeSeries of sensitive distance (range)
from .range import (
    _metadata,
    _segments_diff,
    get_range_channel,
    get_range,
    get_range_spectrogram,
    get_range_spectrum,
)
