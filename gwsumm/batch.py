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

"""Pipeline generator for the Gravitational-wave interferometer
summary information system (`gwsumm`)

This module constructs a directed acyclic graph (DAG) that defines
a workflow to be submitted via the HTCondor scheduler.
"""

import argparse
import os
import shutil
import sys

from glue import pipeline

from gwdatafind.utils import find_credential

from gwpy.io import kerberos as gwkerberos

from gwdetchar import cli

from . import __version__
from .utils import mkdir

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__credits__ = 'Alex Urban <alexander.urban@ligo.org>'

PROG = ('python -m gwsumm.batch' if sys.argv[0].endswith('.py')
        else os.path.basename(sys.argv[0]))


# -- utilities ----------------------------------------------------------------

class GWSummaryJob(pipeline.CondorDAGJob):
    """Job representing a configurable instance of gw_summary.
    """
    logtag = '$(cluster)-$(process)'

    def __init__(self, universe, tag='gw_summary',
                 subdir=None, logdir=None, **cmds):
        pipeline.CondorDAGJob.__init__(self, universe, sys.executable)
        if subdir:
            subdir = os.path.abspath(subdir)
            self.set_sub_file(os.path.join(subdir, '%s.sub' % (tag)))
        if logdir:
            logdir = os.path.abspath(logdir)
            self.set_log_file(os.path.join(
                logdir, '%s-%s.log' % (tag, self.logtag)))
            self.set_stderr_file(os.path.join(
                logdir, '%s-%s.err' % (tag, self.logtag)))
            self.set_stdout_file(os.path.join(
                logdir, '%s-%s.out' % (tag, self.logtag)))
        cmds.setdefault('getenv', 'True')
        for key, val in cmds.items():
            if hasattr(self, 'set_%s' % key.lower()):
                getattr(self, 'set_%s' % key.lower())(val)
            else:
                self.add_condor_cmd(key, val)
        # add python module sub-command
        self._command = ' '.join(['-m', __package__])

    def add_opt(self, opt, value=''):
        pipeline.CondorDAGJob.add_opt(self, opt, str(value))
    add_opt.__doc__ = pipeline.CondorDAGJob.add_opt.__doc__

    def set_command(self, command):
        self._command = ' '.join([
            self._command,
            command,
        ])

    def get_command(self):
        return self._command

    def write_sub_file(self):
        pipeline.CondorDAGJob.write_sub_file(self)
        # insert positional arguments in the right place
        with open(self.get_sub_file(), 'r') as f:
            sub = f.read()
        sub = sub.replace(
            'arguments = "',
            'arguments = " {0}'.format(self.get_command()),
        )
        with open(self.get_sub_file(), 'w') as f:
            f.write(sub)


class GWSummaryDAGNode(pipeline.CondorDAGNode):
    def get_cmd_line(self):
        # merge positional arguments with options
        return ' '.join([
            self.job().get_command(),
            pipeline.CondorDAGNode.get_cmd_line(self),
        ])


# -- parse command-line -------------------------------------------------------

class GWHelpFormatter(argparse.HelpFormatter):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault('indent_increment', 4)
        super(GWHelpFormatter, self).__init__(*args, **kwargs)


def create_parser():
    """Create a command-line parser for this entry point
    """
    # initialize argument parser
    usage = ('%(prog)s --global-config defaults.ini --config-file '
             'myconfig.ini [--config-file myconfig2.ini] [options]')
    parser = argparse.ArgumentParser(
        prog=PROG,
        usage=usage,
        description=__doc__,
        formatter_class=GWHelpFormatter,
    )
    bopts = parser.add_argument_group("Basic options")
    htcopts = parser.add_argument_group("Condor options")
    copts = parser.add_argument_group(
        "Configuration options",
        "Each --global-config file will be used in all nodes of the workflow, "
        "while a single node will be created for each other --config-file",
    )
    popts = parser.add_argument_group(
        "Process options",
        "Configure how this summary will be processed.",
    )
    outopts = parser.add_argument_group("Output options")
    topts = parser.add_argument_group(
        "Time mode options",
        "Choose a stadard time mode, or a GPS [start, stop) interval",
    )

    # general arguments
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help="show verbose output, default: %(default)s",
    )
    parser.add_argument(
        "-V",
        "--version",
        action="version",
        help="show program's version number and exit",
    )
    parser.version = __version__

    # basic options
    bopts.add_argument(
        '-i',
        '--ifo',
        action='store',
        type=str,
        metavar='IFO',
        help="Instrument to process. If this option is set in "
             "the [DEFAULT] of any of the INI files, giving it "
             "here is redundant.",
    )
    wrapgroup = bopts.add_mutually_exclusive_group()
    wrapgroup.add_argument(
        '-w',
        '--skip-html-wrapper',
        action='store_true',
        default=False,
        help="Do not configure first job for HTML htmlnode, default: "
             "%(default)s. Useful for separating large summary pipeline "
             "across multiple DAGs",
    )
    wrapgroup.add_argument(
        '-W',
        '--html-wrapper-only',
        action='store_true',
        help="Only run first job for HTML htmlnode.",
    )
    bopts.add_argument(
        '-t',
        '--file-tag',
        action='store',
        type=str,
        default='gw_summary_pipe',
        help="file tag for pipeline files, default: %(default)s",
    )

    # HTCondor options
    htcopts.add_argument(
        '-u',
        '--universe',
        action='store',
        type=str,
        default='vanilla',
        help="Universe for condor jobs, default: %(default)s",
    )
    htcopts.add_argument(
        '-l',
        '--log-dir',
        action='store',
        type=str,
        default=os.environ.get('LOCALDIR', None),
        help="Directory path for condor log files, default: %(default)s",
    )
    htcopts.add_argument(
        '-m',
        '--maxjobs',
        action='store',
        type=int,
        default=None,
        metavar='N',
        help="Restrict the DAG to submit only N jobs at any one "
             "time, default: %(default)s",
    )
    htcopts.add_argument(
        '-T',
        '--condor-timeout',
        action='store',
        type=float,
        default=None,
        metavar='T',
        help='Configure condor to terminate jobs after T hours '
             'to prevent idling, default: %(default)s',
    )
    htcopts.add_argument(
        '-c',
        '--condor-command',
        action='append',
        type=str,
        default=[],
        help="Extra condor submit commands to add to gw_summary submit file. "
             "Can be given multiple times in the form \"key=value\"",
    )

    # configuration options
    copts.add_argument(
        '-f',
        '--config-file',
        action='append',
        type=str,
        metavar='FILE',
        default=[],
        help="INI file for analysis, may be given multiple times",
    )
    copts.add_argument(
        '-g',
        '--global-config',
        action='append',
        type=str,
        metavar='FILE',
        default=[],
        help="INI file for use in all workflow jobs, may be given "
             "multiple times",
    )
    copts.add_argument(
        '-p',
        '--priority',
        action='append',
        type=str,
        default=[],
        help="priority for DAG node, should be given once "
             "for each --config-file in the same order",
    )

    # process options
    popts.add_argument(
        '--nds',
        action='store_true',
        default='guess',
        help='use NDS as the data source, default: %(default)s',
    )
    popts.add_argument(
        '--single-process',
        action='store_true',
        default=False,
        help="restrict gw_summary to a single process, mainly for "
             "debugging purposes, default: %(default)s",
    )
    popts.add_argument(
        '--multi-process',
        action='store',
        type=int,
        default=None,
        help="maximum number of concurrent sub-processes for each "
             "gw_summary job, {number of CPUs} / {min(number of jobs, 4)}",
    )
    popts.add_argument(
        '-a',
        '--archive',
        action='store_true',
        default=False,
        help="Read archived data from the FILE, and "
             "write back to it at the end",
    )
    popts.add_argument(
        '-S',
        '--on-segdb-error',
        action='store',
        type=str,
        default='raise',
        choices=['raise', 'ignore', 'warn'],
        help="action upon error fetching segments from SegDB",
    )
    popts.add_argument(
        '-G',
        '--on-datafind-error',
        action='store',
        type=str,
        default='raise',
        choices=['raise', 'ignore', 'warn'],
        help="action upon error querying for frames from the "
             "datafind server: default: %(default)s",
    )
    popts.add_argument(
        '--data-cache',
        action='append',
        default=[],
        help='path to LAL-format cache of TimeSeries data files',
    )
    popts.add_argument(
        '--event-cache',
        action='append',
        default=[],
        help='path to LAL-format cache of event trigger files',
    )
    popts.add_argument(
        '--segment-cache',
        action='append',
        default=[],
        help='path to LAL-format cache of state '
             'or data-quality segment files',
    )
    popts.add_argument(
        '--no-htaccess',
        action='store_true',
        default=False,
        help='tell gw_summary to not write .htaccess files',
    )

    # output options
    outopts.add_argument(
        '-o',
        '--output-dir',
        action='store',
        type=str,
        metavar='OUTDIR',
        default=os.curdir,
        help="Output directory for summary information, "
             "default: '%(default)s'",
    )

    # time mode options
    topts.add_argument(
        "--day",
        action="store",
        type=str,
        metavar='YYYYMMDD',
        help="UTC date to process",
    )
    topts.add_argument(
        "--week",
        action="store",
        type=str,
        metavar="YYYYMMDD",
        help="week to process (by UTC starting date)",
    )
    topts.add_argument(
        "--month",
        action="store",
        type=str,
        metavar="YYYYMM",
        help="calendar month to process",
    )
    topts.add_argument(
        "--year",
        action="store",
        type=str,
        metavar="YYYY",
        help="calendar year to process",
    )
    topts.add_argument(
        "-s",
        "--gps-start-time",
        action="store",
        type=int,
        metavar="GPSSTART",
        help="GPS start time",
    )
    topts.add_argument(
        "-e",
        "--gps-end-time",
        action="store",
        type=int,
        metavar="GPSEND",
        help="GPS end time",
    )

    # return the argument parser
    return parser


# -- main code block ----------------------------------------------------------

def main(args=None):
    """Run the command-line Omega scan tool in batch mode
    """
    parser = create_parser()
    args = parser.parse_args(args=args)

    # initialize logger
    logger = cli.logger(
        name=PROG.split('python -m ').pop(),
        level='DEBUG' if args.verbose else 'INFO',
    )

    # check time options
    N = sum([args.day is not None, args.month is not None,
             args.gps_start_time is not None, args.gps_end_time is not None])
    if N > 1 and not (args.gps_start_time and args.gps_end_time):
        raise parser.error("Please give only one of --day, --month, or "
                           "--gps-start-time and --gps-end-time.")

    for (i, cf) in enumerate(args.config_file):
        args.config_file[i] = ','.join(map(os.path.abspath, cf.split(',')))
    args.global_config = list(map(
        os.path.abspath,
        [fp for csv in args.global_config for fp in csv.split(',')],
    ))

    # -- build workflow directories -----------------

    # move to output directory
    indir = os.getcwd()
    mkdir(args.output_dir)
    os.chdir(args.output_dir)
    outdir = os.curdir

    # set node log path, and condor log path
    logdir = os.path.join(outdir, 'logs')
    htclogdir = args.log_dir or logdir
    mkdir(logdir, htclogdir)

    # set config directory and copy config files
    etcdir = os.path.join(outdir, 'etc')
    mkdir(etcdir)

    for (i, fp) in enumerate(args.global_config):
        inicopy = os.path.join(etcdir, os.path.basename(fp))
        if not os.path.isfile(inicopy) or not os.path.samefile(fp, inicopy):
            shutil.copyfile(fp, inicopy)
        args.global_config[i] = os.path.abspath(inicopy)
    for (i, csv) in enumerate(args.config_file):
        inicopy = []
        for fp in csv.split(','):
            fp2 = os.path.join(etcdir, os.path.basename(fp))
            if not os.path.isfile(fp2) or not os.path.samefile(fp, fp2):
                shutil.copyfile(fp, fp2)
            inicopy.append(os.path.abspath(fp2))
        args.config_file[i] = ','.join(inicopy)
    logger.debug("Copied all INI configuration files to %s" % etcdir)

    # -- configure X509 and kerberos for condor -----

    if args.universe != 'local':
        # copy X509 grid certificate into local location
        (x509cert, _) = find_credential()
        x509copy = os.path.join(etcdir, os.path.basename(x509cert))
        shutil.copyfile(x509cert, x509copy)

        # rerun kerberos with new path
        krb5cc = os.path.abspath(os.path.join(etcdir, 'krb5cc.krb5'))
        gwkerberos.kinit(krb5ccname=krb5cc)
        logger.debug("Configured Condor and Kerberos "
                     "for NFS-shared credentials")

    # -- build DAG ----------------------------------

    dag = pipeline.CondorDAG(os.path.join(htclogdir, '%s.log' % args.file_tag))
    dag.set_dag_file(os.path.join(outdir, args.file_tag))

    universe = args.universe

    # -- parse condor commands ----------------------

    # parse into a dict
    condorcmds = {}
    if args.condor_timeout:
        condorcmds['periodic_remove'] = (
            'CurrentTime-EnteredCurrentStatus > %d' %
            (3600 * args.condor_timeout)
        )
    for cmd_ in args.condor_command:
        (key, value) = cmd_.split('=', 1)
        condorcmds[key.rstrip().lower()] = value.strip()

    if args.universe != 'local':
        # add X509 to environment
        for (env_, val_) in zip(['X509_USER_PROXY', 'KRB5CCNAME'],
                                [os.path.abspath(x509copy), krb5cc]):
            condorenv = '%s=%s' % (env_, val_)
            if ('environment' in condorcmds and
                    env_ not in condorcmds['environment']):
                condorcmds['environment'] += ';%s' % condorenv
            elif 'environment' not in condorcmds:
                condorcmds['environment'] = condorenv

    # -- build individual gw_summary jobs -----------

    globalconfig = ','.join(args.global_config)

    jobs = []
    if not args.skip_html_wrapper:
        htmljob = GWSummaryJob(
            'local', subdir=outdir, logdir=logdir,
            tag='%s_local' % args.file_tag, **condorcmds)
        jobs.append(htmljob)
    if not args.html_wrapper_only:
        datajob = GWSummaryJob(
            universe, subdir=outdir, logdir=logdir,
            tag=args.file_tag, **condorcmds)
        jobs.append(datajob)

    # add common command-line options
    for job in jobs:
        if args.day:
            job.set_command('day')
            job.add_arg(args.day)
        elif args.week:
            job.set_command('week')
            job.add_arg(args.week)
        elif args.month:
            job.set_command('month')
            job.add_arg(args.month)
        elif args.year:
            job.set_command('year')
            job.add_arg(args.year)
        elif args.gps_start_time or args.gps_end_time:
            job.set_command('gps')
            job.add_arg(str(args.gps_start_time))
            job.add_arg(str(args.gps_end_time))
        else:
            job.set_command('day')
        if args.nds is True:
            job.add_opt('nds')
        if args.single_process:
            job.add_opt('single-process')
        elif args.multi_process is not None:
            job.add_opt('multi-process', args.multi_process)
        if args.verbose:
            job.add_opt('verbose')
        if args.ifo:
            job.add_opt('ifo', args.ifo)
        job.add_opt('on-segdb-error', args.on_segdb_error)
        job.add_opt('on-datafind-error', args.on_datafind_error)
        job.add_opt('output-dir', outdir)
        for (opt, fplist) in zip(
                ['--data-cache', '--event-cache', '--segment-cache'],
                [args.data_cache, args.event_cache, args.segment_cache]):
            if fplist:
                job.add_arg('%s %s' % (opt, (' %s ' % opt).join(fplist)))
        if args.no_htaccess:
            job.add_opt('no-htaccess')

    # make surrounding HTML first
    if not args.skip_html_wrapper:
        htmljob.add_opt('html-only', '')
        htmljob.add_opt('config-file', ','.join(
            [globalconfig]+args.config_file).strip(','))

        htmlnode = GWSummaryDAGNode(htmljob)
        for configfile in args.config_file:
            htmlnode.add_input_file(args.config_file)
        htmlnode.set_category('gw_summary')
        dag.add_node(htmlnode)
        logger.debug(" -- Configured HTML htmlnode job")

    # create node for each config file
    if not args.html_wrapper_only:
        # add html opts
        datajob.add_opt('no-html', '')
        if args.archive:
            datajob.add_condor_cmd('+SummaryNodeType', '"$(macroarchive)"')
        # configure each data node
        for (i, configfile) in enumerate(args.config_file):
            node = GWSummaryDAGNode(datajob)
            node.add_var_arg('--config-file %s' % ','.join(
                [globalconfig, configfile]).strip(','))
            if args.archive:
                jobtag = os.path.splitext(os.path.basename(configfile))[0]
                archivetag = jobtag.upper().replace('-', '_')
                if args.ifo and archivetag.startswith('%s_' %
                                                      args.ifo.upper()):
                    archivetag = archivetag[3:]
                node.add_var_opt('archive', archivetag)
            for cf in configfile.split(','):
                node.add_input_file(cf)
            node.set_category('gw_summary')
            try:
                node.set_priority(args.priority[i])
            except IndexError:
                node.set_priority(0)
            node.set_retry(1)
            if not args.skip_html_wrapper:
                node.add_parent(htmlnode)
            dag.add_node(node)
            logger.debug(" -- Configured job for config %s" % configfile)

    if args.maxjobs:
        dag.add_maxjobs_category('gw_summary', args.maxjobs)

    # -- finish up ----------------------------------

    dag.write_sub_files()
    dag.write_dag()
    dag.write_script()
    logger.debug("Setup complete, DAG written to: {}".format(
            os.path.abspath(dag.get_dag_file())))

    # return to original directory
    os.chdir(indir)


# -- run from command-line ----------------------------------------------------

if __name__ == "__main__":
    main()
