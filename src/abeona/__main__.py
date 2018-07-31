import argparse
import sys
import logging
from . import __version__

logger = logging.getLogger('abeona')


def main_without_argv():
    argv = sys.argv[1:]
    return main(argv)


def main(argv):
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
    args = parser.parse_args(argv)

    return subcommands[args.sub_command](args.args)


def assemble_main(argv):
    from .assemble.cli import main
    return main(argv)


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

    args = parser.parse_args(args=argv)

    from .subgraphs import main
    return main(args)


if __name__ == '__main__':
    exit(main_without_argv())
