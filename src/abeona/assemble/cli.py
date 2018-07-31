import json
from pathlib import Path
import os


def validate_args(args, parser):
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
    for file_arg in [
        'fastx_forward', 'fastx_reverse', 'fastx_single', 'initial_contigs',
        'kallisto_fastx_forward', 'kallisto_fastx_reverse', 'kallisto_fastx_single'
    ]:
        file = getattr(args, file_arg)
        if file is not None:
            file = Path(file)
            assert file.is_file()
            file = file.absolute()
            assert file.is_file()
            setattr(args, file_arg, str(file))


def main(argv):
    import argparse
    parser = argparse.ArgumentParser(description='Run abeona assembly pipeline.',
                                     prog='abeona assemble')
    parser.add_argument('-o', '--out-dir', help='Output directory', required=True)
    parser.add_argument('-j', '--jobs', help='Number of jobs to schedule concurrently', default=2)
    parser.add_argument('-k', '--kmer-size', default=47)
    parser.add_argument('-m', '--memory', default=3, help='Maximum memory to use in giga bytes')
    parser.add_argument('-q', '--quiet', action='store_true', default=False)

    parser.add_argument('--fastx-forward', help='Forward sequences in FASTA/FASTQ format',
                        required=False)
    parser.add_argument('--fastx-reverse', help='Reverse sequences in FASTA/FASTQ format',
                        required=False)
    parser.add_argument('--fastx-single', help='Single-end sequences in FASTA/FASTQ format',
                        required=False)
    parser.add_argument('--initial-contigs', help='Only start assembly from contigs in this FASTA',
                        required=False)

    group = parser.add_argument_group('graph traversal cleaning arguments')
    group.add_argument('--min-tip-length', type=int, default=0)
    group.add_argument('--min-unitig-coverage', type=int, default=4)

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
    validate_args(args, parser)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(exist_ok=True)

    package_snakefile = Path(__file__).parent / 'Snakefile'
    snakefile = out_dir / 'Snakefile'
    if snakefile.is_file():
        snakefile.unlink()
    snakefile.symlink_to(package_snakefile)

    with open(out_dir / 'args.json', 'w') as fh:
        json.dump({a: getattr(args, a) for a in dir(args) if not a.startswith('_')}, fh)
    cmd = f'cd {out_dir} && snakemake -j {args.jobs}'
    if args.quiet:
        cmd += ' --quiet'
    os.system(cmd)
