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
    from .subgraphs import main
    return main(argv)


if __name__ == '__main__':
    exit(main_without_argv())
