import sys

from .subgraphs import main


def main_without_argv():
    return main(sys.argv[1:])


if __name__ == '__main__':
    exit(main_without_argv())
