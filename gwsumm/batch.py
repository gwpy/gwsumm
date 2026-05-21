# -*- coding: utf-8 -*-
# Copyright (C) Duncan Macleod (2013)
#               Evan Goetz (2026)
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

import os
from pathlib import Path
import shutil
import sys

from glue import pipeline

from gwdetchar.utils import cli
from gwpy.time import tconvert

from . import mode
from .utils import create_parser

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__credits__ = ('Alex Urban <alexander.urban@ligo.org>, '
               'Evan Goetz <evan.goetz@ligo.org>, '
               'Iara Ota <iara.ota@ligo.org>'
               )


PROG = ('python -m gwsumm.batch' if sys.argv[0].endswith('.py')
        else os.path.basename(sys.argv[0]))


# -- utilities ----------------------------------------------------------------

class GWSummaryJob(pipeline.CondorDAGJob):
    """Job representing a configurable instance of gw_summary.
    """
    logtag = '$(cluster)-$(process)'

    def __init__(self, universe, executable='python3', tag='gw_summary',
                 subdir=None, logdir=None, **cmds):
        pipeline.CondorDAGJob.__init__(self, universe, executable)
        if subdir:
            subdir = os.path.abspath(subdir)
            self.set_sub_file(os.path.join(subdir, '%s.sub' % tag))
        if logdir:
            logdir = os.path.abspath(logdir)
            self.set_log_file(os.path.join(
                logdir, '%s-%s.log' % (tag, self.logtag)))
            self.set_stderr_file(os.path.join(
                logdir, '%s-%s.err' % (tag, self.logtag)))
            self.set_stdout_file(os.path.join(
                logdir, '%s-%s.out' % (tag, self.logtag)))
        for key, val in cmds.items():
            if hasattr(self, 'set_%s' % key.lower()):
                getattr(self, 'set_%s' % key.lower())(val)
            else:
                self.add_condor_cmd(key, val)
        # add python module sub-command
        if executable != 'python3':
            self._command = 'exec'
        else:
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


# -- main code block ----------------------------------------------------------

def main(args=None):
    """Run the command-line Omega scan tool in batch mode
    """
    parser = create_parser(include_condor_opts=True)
    args = parser.parse_args(args=args)

    # initialize logger
    logger = cli.logger(
        name=PROG.split('python -m ').pop(),
        level='DEBUG' if args.verbose else 'INFO',
    )

    # set mode and output directory
    # This will set the output directory to be day/YYYYMMDD, week/YYYYMMDD,
    # month/YYYYMM, or gps/<gps start>-<gps end>.
    mode.set_mode(args.mode)
    try:
        utc = tconvert(args.gpsstart)
        outpath = mode.get_base(utc)
        outpath = str(Path(outpath).parent)
    except ValueError:
        outpath = 'gps'

    for (i, cf) in enumerate(args.config_file):
        args.config_file[i] = ','.join(map(os.path.abspath, cf.split(',')))
    args.global_config = list(map(
        os.path.abspath,
        [fp for csv in args.global_config for fp in csv.split(',')],
    ))

    # -- build workflow directories -----------------

    # move to output directory
    indir = os.getcwd()
    os.makedirs(args.output_dir, exist_ok=True)
    os.chdir(args.output_dir)
    outdir = os.curdir

    # set node log path, and condor log path
    logdir = os.path.join(outdir, 'logs')
    htclogdir = args.log_dir or logdir
    os.makedirs(logdir, exist_ok=True)
    os.makedirs(htclogdir, exist_ok=True)

    # set config directory and copy config files
    etcdir = os.path.join(outdir, 'etc')
    os.makedirs(etcdir, exist_ok=True)

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

    # -- build DAG ----------------------------------

    dag = pipeline.CondorDAG(os.path.join(htclogdir, f'{args.file_tag}.log'))
    dag.set_dag_file(os.path.join(outdir, args.file_tag))

    # -- parse condor commands ----------------------

    # parse into a dict
    condorcmds = {}
    if args.condor_timeout:
        condorcmds['periodic_remove'] = (
            f'CurrentTime-EnteredCurrentStatus > {3600 * args.condor_timeout}'
        )
    for cmd_ in args.condor_command:
        (key, value) = cmd_.split('=', 1)
        # TODO: the next lines either add the key-value to the condorcmds
        #  dictionary or they *append* the value to the existing value already
        #  in the condorcmds dictionary. It's not a preferred solution, but
        #  the reason it is done this way is because of the way the
        #  LIGO summary pages configuration scripts work (the echo_and_run
        #  command). We would need to decide if the scripts or function
        #  should be rewritten
        if key in condorcmds:
            logger.debug(f"Appending {value.strip()} to condor '{key}' value")
            condorcmds[key] += f" {value.strip()}"
        else:
            condorcmds[key.rstrip().lower()] = value.strip()

    # Use scitokens
    condorcmds['use_oauth_services'] = 'scitokens'
    if ('environment' in condorcmds and
            'BEARER_TOKEN_FILE' not in condorcmds['environment']):
        condorcmds['environment'] += (
            ' BEARER_TOKEN_FILE='
            '$$(CondorScratchDir)/.condor_creds/scitokens.use'
            )
    elif 'environment' not in condorcmds:
        condorcmds['environment'] = (
            'BEARER_TOKEN_FILE='
            '$$(CondorScratchDir)/.condor_creds/scitokens.use'
            )

    # set XDG_CONFIG_HOME and XDG_CACHE_HOME vars for astropy
    # see https://docs.astropy.org/en/stable/environment_variables.html
    # the logic here is that if 'environment' exists in condorcmds then
    # if XDG_CONFIG_HOME or XDG_CACHE_HOME doesn't exist then add it (space
    # separated list), or if 'environment' doesn't exist then we've already
    # started with BEARER_TOKEN_FILE above and now we add the XDG variables
    if (('environment' in condorcmds and
         'XDG_CONFIG_HOME' not in condorcmds['environment']) or
            'environment' not in condorcmds):
        condorcmds['environment'] += (
            ' XDG_CONFIG_HOME=$$(CondorScratchDir)/tmp/config'
            )
    if (('environment' in condorcmds and
         'XDG_CACHE_HOME' not in condorcmds['environment']) or
            'environment' not in condorcmds):
        condorcmds['environment'] += (
            ' XDG_CACHE_HOME=$$(CondorScratchDir)/tmp/cache'
            )

    # Make sure the environment is a quotation marked string
    if not condorcmds['environment'].startswith('"'):
        condorcmds['environment'] = f"\"{condorcmds['environment']}"
    if not condorcmds['environment'].endswith('"'):
        condorcmds['environment'] = f"{condorcmds['environment']}\""

    # Environment variables from the access point
    envvars = ('DEFAULT_SEGMENT_SERVER, GWDATAFIND_SERVER, NDSSERVER, '
               'CONDA_EXE')

    # -- build individual gw_summary jobs -----------

    globalconfig = ','.join(args.global_config)

    jobs = []
    # Define an HTML only job which runs in the local universe.
    # This is always just one job, based on the structure we've defined.
    if args.html_only or not args.no_html:
        htmljob = GWSummaryJob(
            'local', executable='apptainer',
            subdir=outdir, logdir=logdir,
            tag=f'{args.file_tag}_local', **condorcmds,
            getenv=envvars,
        )
        jobs.append(htmljob)
    # Define data jobs which run either in the local or container universe. We
    # need to handle the case where jobs may not have completed previously
    # so there may be no archive file.
    if not args.html_only or args.no_html:
        # HTCondor file transfer commands
        transfer_aux_files = []
        if args.universe != 'local':
            # set executable
            executable = 'python3'
            # Set transfer commands
            for cmd, val in {'should_transfer_files': 'YES',
                             'transfer_input_files': '$(inputfiles)',
                             'transfer_output_files': outpath,
                             'container_image': args.container_path}.items():
                if cmd not in condorcmds:
                    condorcmds[cmd] = val
            # Change the environment variable definitions when transferring
            # files to an execute point.
            # List the all of the environment variables set by the submit file
            condor_envvars = condorcmds['environment'].split(' ')
            # When an environment variable points to a directory or file,
            # change the environment variable just to the name and add the
            # path to the list of files to be transferred.
            for idx, entry in enumerate(condor_envvars):
                var, val = entry.split('=', 1)
                if (testpath := Path(val)).exists():
                    transfer_aux_files.append(str(testpath))
                    updated_env = testpath.name
                    condor_envvars[idx] = (f'{var}=$$(CondorScratchDir)/'
                                           f'{updated_env}')
            # Now rejoin the environment variables
            condorcmds['environment'] = ' '.join(condor_envvars)
            # Join the auxiliary files in a comma separated list
            transfer_aux_files = ','.join(transfer_aux_files)

            # container universe jobs condor commands
            for cmd_ in args.condor_command_container:
                (key, value) = cmd_.split('=', 1)
                if key in condorcmds:
                    logger.debug(f"Appending {value.strip()} to condor '{key}'"
                                 f" value")
                    condorcmds[key] += f" {value.strip()}"
                else:
                    condorcmds[key.rstrip()] = value.strip()
        else:
            executable = 'apptainer'
        # A data job
        datajob = GWSummaryJob(
            args.universe, executable=executable,
            subdir=outdir, logdir=logdir,
            tag=args.file_tag, **condorcmds,
            getenv=envvars,
        )
        # A data job with no archive file
        datajob_noarchive = GWSummaryJob(
            args.universe, executable=executable,
            subdir=outdir, logdir=logdir,
            tag=f'{args.file_tag}_noarchive', **condorcmds,
            getenv=envvars,
        )
        jobs.append(datajob)
        jobs.append(datajob_noarchive)

    # add common command-line options
    for job in jobs:
        if job.get_universe() == 'local':
            job.set_command(
                f'{os.path.abspath(args.container_path)} python3 -m '
                f'{__package__}'
            )
        job.set_command(args.mode)
        if args.mode == 'day':
            job.add_arg(args.day.strftime("%Y%m%d"))
        elif args.mode == 'week':
            job.add_arg(args.week.strftime("%Y%m%d"))
        elif args.mode == 'month':
            job.add_arg(args.month.strftime("%Y%m"))
        elif args.mode == 'gps':
            job.add_arg(str(args.gpsstart))
            job.add_arg(str(args.gpsend))
        if args.nds is True:
            job.add_opt('nds')
        if args.multiprocess is not None:
            job.add_opt('multi-process', args.multiprocess)
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
        # If the job is in the container universe, the transferred archive file
        # will be in the flat directory so set the archive read directory for
        # the job to read from the flat directory
        if ('noarchive' not in job.get_sub_file() and
                args.archive and
                job.get_universe() == 'container'):
            job.add_opt('archive-read-dir', '.')
        if (args.archive and job.get_universe() == 'container' and
                args.archive_write_dir):
            job.add_opt('archive-write-dir', args.archive_write_dir)

    # make surrounding HTML first
    if args.html_only or not args.no_html:
        htmljob.add_opt('html-only', '')
        htmljob.add_opt('config-file', ','.join(
            [globalconfig]+args.config_file).strip(','))

        htmlnode = GWSummaryDAGNode(htmljob)
        htmlnode.add_input_file(args.config_file)
        htmlnode.set_category('gw_summary')
        dag.add_node(htmlnode)
        logger.debug(" -- Configured HTML htmlnode job")

    # create node for each config file
    if not args.html_only or args.no_html:
        # add html opts
        datajob.add_opt('no-html', '')
        datajob_noarchive.add_opt('no-html', '')
        if args.archive:
            datajob.add_condor_cmd('+SummaryNodeType', '"$(macroarchive)"')
            datajob_noarchive.add_condor_cmd(
                '+SummaryNodeType', '"$(macroarchive)"'
            )
        # configure each data node
        for (i, configfile) in enumerate(args.config_file):
            # If we are using archive files, then we need to transfer any
            # HDF5 archive files with the tag in args.archive if running in
            # the container universe
            files = ''  # Comma separated list of existing archive files
            if args.archive:
                # First figure out the archive tag
                jobtag = os.path.splitext(os.path.basename(configfile))[0]
                archivetag = jobtag.upper().replace('-', '_')
                if args.ifo and archivetag.startswith(f'{args.ifo.upper()}_'):
                    archivetag = archivetag[3:]

                # If running in the container universe, check to see if any
                # archive file already exists. If it does then we need to
                # transfer it there.
                # Both datajob types have the same universe, so only check one
                if datajob.get_universe() == 'container':
                    # If we find any matching archive files make a comma
                    # separated string with those files
                    if len(archives := sorted(Path(args.archive_read_dir).glob(
                            f'*-{archivetag}-*.h5'))) > 0:
                        archives = [str(f) for f in archives]
                        files = ','.join(archives).strip(',')

                # If there was no archive file found for this config group and
                # running in the container universe for the data jobs, then
                # this node will be a "no archive". We don't transfer an
                # archive file, and we don't give the --archive-read-dir .
                # option.
                if len(files) == 0 and datajob.get_universe() == 'container':
                    node = GWSummaryDAGNode(datajob_noarchive)
                else:
                    node = GWSummaryDAGNode(datajob)

                # Finally, we have a node, so add the --archive <archivetag>
                # option
                node.add_var_opt('archive', archivetag)
            else:
                node = GWSummaryDAGNode(datajob)

            # Configuration files for gw_summary jobs are going to be
            # given as the full path (local universe jobs), or just the
            # file name (container universe jobs)
            if datajob.get_universe() == 'local':
                config_files = ','.join([globalconfig, configfile]).strip(',')
            elif datajob.get_universe() == 'container':
                # Use just the file name
                config_files = ','.join(
                    [Path(f).name for f in
                     ','.join(
                         [globalconfig.strip(','), configfile.strip(',')]
                     ).split(',')]
                )
                # Transferred files are the paths to file relative to the
                # submit location (usually ~detchar/public_html/summary).
                # The comma separated string is added to the dag as a macro
                # variable inputfiles
                txfr_files = ','.join(
                    [globalconfig, configfile, transfer_aux_files, files]
                ).strip(',')
                node.add_macro('inputfiles', txfr_files)
            else:
                raise ValueError(f"Unknown universe {node.get_universe()}")
            node.add_var_opt('config-file', config_files)
            node.set_category('gw_summary')
            try:
                node.set_priority(args.priority[i])
            except IndexError:
                node.set_priority(0)
            node.set_retry(1)
            if not args.no_html:
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
