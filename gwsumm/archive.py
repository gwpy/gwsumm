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

"""This module handles HDF archiving of data.

In production for LIGO, the LIGO Summary Pages (LSP) Service runs at regular
intervals (about every 10 minutes), so an HDF5 file is used to archive the
data read and produced from one instance so that the next instance doesn't
have to re-read and re-produce the same data.

All data products are stored just using the 'standard' gwpy `.write()` method
for that object.
"""

import tempfile
import shutil
import warnings
import re
import datetime
import os

from numpy import (unicode_, ndarray)

from astropy.table import Table

from gwpy.time import (from_gps, to_gps)
from gwpy.timeseries import (StateVector, TimeSeries)
from gwpy.spectrogram import Spectrogram
from gwpy.segments import (SegmentList, Segment, DataQualityFlag)

from . import (globalv, mode)
from .data import (get_channel, add_timeseries, add_spectrogram,
                   add_coherence_component_spectrogram)
from .triggers import (EventTable, add_triggers)

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'

re_rate = re.compile('_EVENT_RATE_')


def write_data_archive(outfile, channels=True, timeseries=True,
                       spectrogram=True, segments=True, triggers=True):
    """Build and save an HDF archive of data processed in this job.

    Parameters
    ----------
    outfile : `str`
        path to target HDF5 file

    timeseries : `bool`, optional
        include `TimeSeries` data in archive

    spectrogram : `bool`, optional
        include `Spectrogram` data in archive

    segments : `bool`, optional
        include `DataQualityFlag` data in archive

    triggers : `bool`, optional
        include `EventTable` data in archive
    """
    from h5py import File

    # backup existing file, in case something fails
    backup = backup_existing_archive(outfile)

    try:
        with File(outfile, 'w') as h5file:

            # -- channels -----------------------

            if channels and globalv.CHANNELS:
                cols = ('name', 'sample_rate', 'frametype', 'unit')
                rows = []
                for chan in globalv.CHANNELS:
                    rows.append((
                        chan.ndsname,
                        chan.sample_rate.to('Hz').value if
                        chan.sample_rate is not None else 0,
                        str(getattr(chan, 'frametype', None)) or '',
                        str(chan.unit) if chan.unit else '',
                    ))
                Table(names=cols, rows=rows).write(h5file, 'channels')

            # -- timeseries ---------------------

            if timeseries:
                tgroup = h5file.create_group('timeseries')
                sgroup = h5file.create_group('statevector')
                # loop over channels
                for c, tslist in globalv.DATA.items():
                    c = get_channel(c)
                    # ignore trigger rate TimeSeries
                    if re_rate.search(str(c)):
                        continue
                    # loop over time-series
                    for ts in tslist:
                        # ignore fast channels who weren't used
                        # for a timeseries:
                        if (not isinstance(ts, StateVector) and
                                ts.sample_rate.value > 16.01 and
                                not getattr(c, '_timeseries', False)):
                            continue
                        # archive timeseries
                        try:
                            name = '%s,%s,%s' % (ts.name, ts.channel.ndsname,
                                                 ts.t0.value)
                        except AttributeError:
                            name = '%s,%s' % (ts.name, ts.t0.value)
                        if isinstance(ts, StateVector):
                            group = sgroup
                        else:
                            group = tgroup
                        _write_object(ts, group, path=name, format='hdf5')

            # -- spectrogram --------------------

            if spectrogram:
                for tag, gdict in zip(
                        ['spectrogram', 'coherence-components'],
                        [globalv.SPECTROGRAMS, globalv.COHERENCE_COMPONENTS]):
                    group = h5file.create_group(tag)
                    # loop over channels
                    for key, speclist in gdict.items():
                        # loop over time-series
                        for spec in speclist:
                            name = '%s,%s' % (key, spec.t0.value)
                            _write_object(spec, group, path=name,
                                          format='hdf5')

            # -- segments -----------------------

            if segments:
                group = h5file.create_group('segments')
                # loop over channels
                for name, dqflag in globalv.SEGMENTS.items():
                    dqflag.write(group, path=name, format='hdf5')

            # -- triggers -----------------------

            if triggers:
                group = h5file.create_group('triggers')
                for key in globalv.TRIGGERS:
                    archive_table(globalv.TRIGGERS[key], key, group)

    except Exception:  # if it fails for any reason, reinstate the backup
        if backup:
            restore_backup(backup, outfile)
        raise
    finally:  # remove the backup regardless of what happens
        if backup is not None and os.path.isfile(backup):
            os.remove(backup)


def read_data_archive(sourcefile):
    """Read archived data from an HDF5 archive source

    This method reads all found data into the data containers defined by
    the `gwsumm.globalv` module, then returns nothing.

    Parameters
    ----------
    sourcefile : `str`
        path to source HDF5 file
    """
    from h5py import File

    with File(sourcefile, 'r') as h5file:

        # -- channels ---------------------------

        try:
            ctable = Table.read(h5file['channels'])
        except KeyError:  # no channels table written
            pass
        else:
            for row in ctable:
                chan = get_channel(row['name'])
                for p in ctable.colnames[1:]:
                    if row[p]:
                        setattr(chan, p, row[p])

        # -- timeseries -------------------------

        for dataset in h5file.get('timeseries', {}).values():
            ts = TimeSeries.read(dataset, format='hdf5')
            if (re.search(r'\.(rms|min|mean|max|n)\Z', ts.channel.name) and
                    ts.sample_rate.value == 1.0):
                ts.channel.type = 's-trend'
            elif re.search(r'\.(rms|min|mean|max|n)\Z', ts.channel.name):
                ts.channel.type = 'm-trend'
            ts.channel = get_channel(ts.channel)
            try:
                add_timeseries(ts, key=ts.channel.ndsname)
            except ValueError:
                if mode.get_mode() != mode.Mode.day:
                    raise
                warnings.warn('Caught ValueError in combining daily archives')
                # get end time
                globalv.DATA[ts.channel.ndsname].pop(-1)
                t = globalv.DATA[ts.channel.ndsname][-1].span[-1]
                add_timeseries(ts.crop(start=t), key=ts.channel.ndsname)

        # -- statevector -- ---------------------

        for dataset in h5file.get('statevector', {}).values():
            sv = StateVector.read(dataset, format='hdf5')
            sv.channel = get_channel(sv.channel)
            add_timeseries(sv, key=sv.channel.ndsname)

        # -- spectrogram ------------------------

        for tag, add_ in zip(
                ['spectrogram', 'coherence-components'],
                [add_spectrogram, add_coherence_component_spectrogram]):
            for key, dataset in h5file.get(tag, {}).items():
                key = key.rsplit(',', 1)[0]
                spec = Spectrogram.read(dataset, format='hdf5')
                spec.channel = get_channel(spec.channel)
                add_(spec, key=key)

        # -- segments ---------------------------

        for name, dataset in h5file.get('segments', {}).items():
            dqflag = DataQualityFlag.read(h5file, path=dataset.name,
                                          format='hdf5')
            globalv.SEGMENTS += {name: dqflag}

        # -- triggers ---------------------------

        for dataset in h5file.get('triggers', {}).values():
            load_table(dataset)


def backup_existing_archive(filename, suffix='.h5',
                            prefix='gw_summary_archive_', dir=None):
    """Create a copy of an existing archive.
    """
    backup = tempfile.mktemp(suffix=suffix, prefix=prefix, dir=dir)
    try:
        shutil.move(filename, backup)
    except IOError:
        return None
    else:
        return backup


def restore_backup(backup, target):
    """Reinstate a backup copy of the archive.
    """
    shutil.move(backup, target)


def find_daily_archives(start, end, ifo, tag, basedir=os.curdir):
    """Find the daily archives spanning the given GPS [start, end) interval
    """
    archives = []
    s = from_gps(to_gps(start))
    e = from_gps(to_gps(end))
    while s < e:
        daybase = mode.get_base(s, mode=mode.Mode.day)
        ds = to_gps(s)
        s += datetime.timedelta(days=1)
        de = to_gps(s)
        archivedir = os.path.join(basedir, daybase, 'archive')
        arch = os.path.join(archivedir, '%s-%s-%d-%d.h5'
                            % (ifo, tag, ds, de-ds))
        if os.path.isfile(arch):
            archives.append(arch)
    return archives


# -- utility methods --------------------------------------------------------

def _write_object(data, *args, **kwargs):
    """Internal method to write something to HDF5 with error handling
    """
    try:
        return data.write(*args, **kwargs)
    except ValueError as e:
        warnings.warn(str(e))
    except RuntimeError as e:
        if 'name already exists' in str(e).lower():
            warnings.warn(str(e))
        else:
            raise


def segments_to_array(segmentlist):
    """Convert a `SegmentList` to a 2-dimensional `numpy.ndarray`
    """
    out = ndarray((len(segmentlist), 2), dtype=float)
    for i, seg in enumerate(segmentlist):
        out[i] = seg
    return out


def segments_from_array(array):
    """Convert a 2-dimensional `numpy.ndarray` to a `SegmentList`
    """
    out = SegmentList()
    for row in array:
        out.append(Segment(*row))
    return out


def archive_table(table, key, parent):
    """Add a table to the given HDF5 group

    .. warning::

       If the input ``table`` is empty, it will not be archived

    Parameters
    ----------
    table : `~astropy.table.Table`
        the data to archive

    key : `str`
        the path (relative to ``parent``) at which to store the table

    parent : `h5py.Group`
        the h5py group in which to add this dataset

    """
    if len(table) == 0:
        warnings.warn("%r table is empty and will not be archived" % key)
        return
    table.meta.pop('psd', None)  # pycbc_live
    table.meta.pop('loudest', None)  # pycbc_live
    try:
        table.meta['segments'] = segments_to_array(table.meta['segments'])
    except KeyError:
        pass
    for col in table.columns:
        if table[col].dtype.type is unicode_:
            table.replace_column(col, table[col].astype(str))
    table.write(parent, path=key, format='hdf5')
    return key


def load_table(dataset):
    """Read table from the given HDF5 group

    The `EventTable` is read, stored in the memory archive, then returned

    Parameters
    ----------
    dataset : `h5py.Dataset`
        n-dimensional table to load from hdf5

    Returns
    -------
    table : `~gwpy.table.EventTable`
        the table of events loaded from hdf5
    """
    table = EventTable.read(dataset, format='hdf5')
    try:
        table.meta['segments'] = segments_from_array(table.meta['segments'])
    except KeyError:
        table.meta['segments'] = SegmentList()
    add_triggers(table, dataset.name.split('/')[-1])
    return table
