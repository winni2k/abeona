"""
abeona subgraphs

Partition cortex graph into its subgraphs.

Usage: cortexpy [--help] <graph> <out_dir>

Options:
   -h, --help  Display this help message
   <graph>     Cortex graph
   <out_dir>   Output dir of subgraphs
"""
import sys
from pathlib import Path
from docopt import docopt

VERSION_STRING='0.0.1'

def main(argv):
    args = docopt(__doc__, argv=sys.argv, version=VERSION_STRING)

    out_dir = Path(args['<out_dir>'])
    out_dir.mkdir(exist_ok=True)

    input_graph = Path(args['<graph>'])
    if not input_graph.is_file():
        raise Exception(f'Input cortex graph ({input_graph}) does not exist')



if __name__ == '__main__':
    try:
        main(sys.argv)
    except e:
        print(e)
        exit(1)
    exit(0)
