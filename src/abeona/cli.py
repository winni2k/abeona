import logging
import subprocess
import sys

logger = logging.getLogger('abeona')


def main(argv=sys.argv):
    from . import __version__
    import argparse
    subcommands = {
        'assemble': assemble_main,
        'subgraphs': subgraphs_main,
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

    group = parser.add_argument_group('graph traversal cleaning')
    group.add_argument('--min-tip-length', type=int, default=0)
    group.add_argument('--min-unitig-coverage', type=int, default=4)

    group = parser.add_argument_group('candidate transcript creation')
    group.add_argument('--max-paths-per-subgraph', type=int, default=0)
    group.add_argument('--merge-candidates-before-kallisto', action='store_true')

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

    args = parser.parse_args(args=argv)
    make_file_paths_absolute(args)
    validate_assemble_args(args, parser)

    from pathlib import Path
    import json

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
