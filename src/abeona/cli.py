import logging
import sys

logger = logging.getLogger('abeona')


def main(argv=sys.argv):
    from . import __version__
    import argparse
    subcommands = {
        'assemble': assemble_main,
        'subgraphs': subgraphs_main,
        'reads': reads_main,
    }
    parser = argparse.ArgumentParser(description=f'abeona version {__version__}', prog='abeona')
    parser.add_argument('-v', '--version', action='version',
                        version=f'%(prog)s version {__version__}')
    parser.add_argument('sub_command', choices=sorted(subcommands.keys()),
                        help='abeona sub-command')
    parser.add_argument('args', nargs=argparse.REMAINDER, help='sub-command arguments')
    args = parser.parse_args(argv[1:])

    return subcommands[args.sub_command](args.args)


def assemble_main(argv):
    import argparse
    parser = argparse.ArgumentParser(description='Run abeona assembly pipeline.',
                                     prog='abeona assemble',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--out-dir', help='Output directory', required=True)
    parser.add_argument('-j', '--jobs', help='Number of jobs to schedule concurrently', default=2)
    parser.add_argument('-k', '--kmer-size', default=47)
    parser.add_argument('-m', '--memory', default=3, help='Maximum memory to use in giga bytes')
    parser.add_argument('-q', '--quiet', action='store_true')
    parser.add_argument('--resume', action='store_true')

    parser.add_argument('--fastx-forward', help='Forward sequences in FASTA/FASTQ format',
                        required=False)
    parser.add_argument('--fastx-reverse', help='Reverse sequences in FASTA/FASTQ format',
                        required=False)
    parser.add_argument('--fastx-single', help='Single-end sequences in FASTA/FASTQ format',
                        required=False)
    parser.add_argument('--initial-contigs', help='Only start assembly from contigs in this FASTA',
                        required=False)
    parser.add_argument('--extra-start-kmer',
                        help='Disconnect this k-mer from incoming k-mers before '
                             'candidate transcript creation. '
                             'This may be useful when assembling circular genomes. '
                             'Best used with --initial-contigs to make sure the k-mer exists '
                             'in the consistent cortexpy graph.',
                        required=False)

    group = parser.add_argument_group('Graph traversal cleaning')
    group.add_argument('--min-tip-length', type=int, default=0,
                       help='Prune tips shorter than this value')
    group.add_argument('--min-unitig-coverage', type=int, default=4,
                       help="Prune unitigs with mean coverage below this value")
    group.add_argument('--prune-tips-with-mccortex', action='store_true',
                       help="Use Mccortex instead of cortexpy to prune unitigs.")

    group = parser.add_argument_group('Candidate transcript creation')
    group.add_argument('--max-paths-per-subgraph', type=int, default=0,
                       help='Ignore graphs that have more than this number of paths')

    group = parser.add_argument_group('kallisto arguments',
                                      description='Arguments passed directly on to kallisto')
    group.add_argument('--bootstrap-samples', type=int, default=100,
                       help='Number of kallisto bootstrap samples')
    group.add_argument('--kallisto-fastx-forward',
                       help='Forward sequences in FASTA/FASTQ format to use for quantification')
    group.add_argument('--kallisto-fastx-reverse',
                       help='Reverse sequences in FASTA/FASTQ format to use for quantification')
    group.add_argument('--kallisto-fastx-single',
                       help='Single-end sequences in FASTA/FASTQ format to use for quantification')
    group.add_argument('--kallisto-fragment-length', type=float)
    group.add_argument('--kallisto-sd', type=float)
    group.add_argument('--kallisto-threads', type=int, default=2,
                       help='Number of logical cores to assign to a single kallisto quant job. '
                            'Needs to be less than or equal to --jobs')

    group = parser.add_argument_group('Candidate transcript filtering')
    group.add_argument('--estimated-count-threshold', type=float, default=1,
                       help='Threshold over which the estimated transcript count from kallisto is '
                            'counted towards keeping a transcript')
    group.add_argument('--bootstrap-proportion-threshold', type=float, default=0.95,
                       help="Proportion of bootstrap iterations for which a transcript's estimated "
                            "counts must be above the --estimated-count-threshold")

    args = parser.parse_args(args=argv)
    make_file_paths_absolute(args)
    validate_assemble_args(args, parser)

    from pathlib import Path
    import json
    import subprocess

    out_dir = Path(args.out_dir)
    out_dir.mkdir(exist_ok=True)

    script_name = 'assemble.nf'
    pipeline_script = out_dir / script_name
    package_script = Path(__file__).parent / 'assemble' / script_name
    if pipeline_script.is_symlink():
        pipeline_script.unlink()
    pipeline_script.symlink_to(package_script)

    args_file = 'args.json'
    args_dict = {a: getattr(args, a) for a in dir(args) if not a.startswith('_')}
    args_dict['mccortex'] = f'mccortex {args.kmer_size}'
    args_dict['mccortex_args'] = f'--sort --force -m {args.memory}G'
    with open(out_dir / args_file, 'w') as fh:
        json.dump(args_dict, fh)
    cmd = f'cd {out_dir} && nextflow run {script_name} -process.maxForks {args.jobs} -params-file {args_file}'
    if args.resume:
        cmd += ' -resume'
    logger.info(cmd)
    return subprocess.run(cmd, shell=True).returncode


def subgraphs_main(argv):
    import argparse
    parser = argparse.ArgumentParser(description='Partition cortex graph into its subgraphs.',
                                     prog='abeona subgraphs')
    parser.add_argument('graph')
    parser.add_argument('out_dir')
    parser.add_argument('-m', '--memory', type=int, default=3)
    parser.add_argument('-c', '--cores', type=int, default=2)
    parser.add_argument('--initial-contigs', help='Only start assembly from contigs in this FASTA',
                        required=False)

    args = parser.parse_args(argv)

    from .subgraphs import main
    return main(args)


def reads_main(argv):
    import argparse
    parser = argparse.ArgumentParser(description='Assign reads to cortex graphs.',
                                     prog='abeona reads')
    parser.add_argument('graph_list',
                        help="""
                        Tab delimited text file with two columns and optional header line.
                        All lines starting with '#' are ignored. For example:
                        ```
                        # prefix\tgraph
                        out_dir/g1\tg1.ctx
                        out_dir/g2\tg2.ctx
                        ```
                        """)
    parser.add_argument('forward', help='Forward reads file in FASTA or FASTQ format.')
    parser.add_argument('--reverse',
                        help='Reverse reads file in FASTA or FASTQ format. '
                             'Only specified if reads are paired',
                        default=None)
    parser.add_argument('--format', help='File format of input reads.',
                        choices=['fasta', 'fastq'],
                        required=True)

    args = parser.parse_args(argv)

    from .reads import main
    return main(args)


def validate_assemble_args(args, parser):
    if (args.fastx_forward is None) != (args.fastx_reverse is None):
        raise parser.error(
            'Both --fastx-forward and --fastx-reverse must be specified or not specified')
    if not args.fastx_forward and not args.fastx_single:
        raise parser.error('Need to specify --fastx-forward and/or --fastx-single')

    if (args.kallisto_fastx_forward is None) != (args.kallisto_fastx_reverse is None):
        raise parser.error('Both --kallisto-fastx-forward and --kallisto-fastx-reverse'
                           ' must be specified or not specified')
    if (args.kallisto_fastx_forward is None) == (args.kallisto_fastx_single is None):
        raise parser.error('Need to specify only --kallisto-fastx-forward'
                           ' or --kallisto-fastx-single')
    if args.kallisto_fastx_single is not None:
        if args.kallisto_fragment_length is None:
            raise parser.error('Required: --kallisto-fragment-length')
        if args.kallisto_sd is None:
            raise parser.error('Required: --kallisto-sd')
    if int(args.jobs) < int(args.kallisto_threads):
        raise parser.error('--jobs needs to be greater or equal to --kallisto-threads')


def make_file_paths_absolute(args):
    from pathlib import Path
    for file_arg in [
        'fastx_forward', 'fastx_reverse', 'fastx_single', 'initial_contigs',
        'kallisto_fastx_forward', 'kallisto_fastx_reverse', 'kallisto_fastx_single'
    ]:
        file = getattr(args, file_arg)
        if file is not None:
            file = Path(file)
            assert file.is_file(), f'Not a file: {file}'
            file = file.absolute()
            assert file.is_file(), f'Not a file: {file}'
            setattr(args, file_arg, str(file))
