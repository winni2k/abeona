"""
abeona subgraphs

Partition cortex graph into its subgraphs.

Usage: abeona subgraphs [options] <graph> <out_dir>

Options:
   -h, --help  Display this help message
   <graph>     Cortex graph
   <out_dir>   Output dir of subgraphs

   -m, --memory <n>  Max memory to use for mccortex (in gigabytes). [default: 3]
   -c, --cores <n>   Max number of cores to use. [default: 2]
"""
import attr
import collections
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from cortexpy.graph.parser.streaming import kmer_list_generator_from_stream
import logging
from subprocess import check_call

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('abeona.subgraphs')


@attr.s(slots=True)
class MccortexRunner(object):
    dist = attr.ib()
    mem = attr.ib()
    mccortex_bin = attr.ib('mccortex')
    threads = attr.ib(2)
    sort = attr.ib(True)

    def build_subgraph(self, *, input, out, initial_kmer, graph_id):
        import shutil

        seed_file = out.with_suffix('.seed.fa')
        SeqIO.write([SeqRecord(Seq(initial_kmer), id=f'{graph_id}_seed', description='')],
                    str(seed_file), 'fasta')
        tmp_out = str(out) + '.tmp'
        kmer_size = str(len(initial_kmer))
        mem_str = f'{self.mem}G'
        check_call(
            ['mccortex', kmer_size, 'subgraph',
             '--force',
             '--out', tmp_out,
             '-m', mem_str,
             '--threads', str(self.threads),
             '--dist', str(self.dist),
             '--seq', str(seed_file),
             '--quiet',
             str(input)]
        )
        check_call(
            ['mccortex', kmer_size, 'sort', '--force', '--quiet',
             '-m', mem_str,
             str(tmp_out)]
        )
        shutil.move(tmp_out, out)


def remove_subgraph_kmers_from_dict(graph, dict):
    """Reads kmer strings from graph and removes the strings from dict in place"""
    with open(graph, 'rb') as fh:
        for kmer_list in kmer_list_generator_from_stream(fh):
            del dict[''.join(kmer_list)]


def main(argv):
    import argparse
    parser = argparse.ArgumentParser(description='Partition cortex graph into its subgraphs.', prog='abeona subgraphs')
    parser.add_argument('graph')
    parser.add_argument('out_dir')
    parser.add_argument('-m', '--memory', type=int, default=3)
    parser.add_argument('-c', '--cores', type=int, default=2)
    parser.add_argument('--initial-contigs', help='Only start assembly from contigs in this FASTA',
                        required=False)

    args = parser.parse_args(args=argv)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(exist_ok=True)

    input_graph = Path(args.graph)
    if not input_graph.is_file():
        raise Exception(f'Input cortex graph ({input_graph}) does not exist')

    logger.info('Loading graph kmer strings')
    unseen_kmer_strings = collections.OrderedDict()
    with open(input_graph, 'rb') as fh:
        for kmer_list in kmer_list_generator_from_stream(fh):
            unseen_kmer_strings[''.join(kmer_list)] = None

    if args.initial_contigs:
        from cortexpy.utils import lexlo
        logger.info('Subsetting graph kmer strings on initial-contig kmers')
        initial_kmers = set()
        for rec in SeqIO.parse(args.initial_contigs,'fasta'):
            kmer_size = len(next(iter(unseen_kmer_strings)))
            assert len(rec.seq) >= kmer_size
            for start_idx in range(len(rec.seq) - kmer_size + 1):
                initial_kmers.add(str(lexlo(rec.seq[start_idx:start_idx+kmer_size])))
        logger.debug(f'initial kmers: {initial_kmers}')

        logger.info(f'Keeping {len(initial_kmers)} out of {len(unseen_kmer_strings)} kmers')
        old_unseen_kmer_strings = unseen_kmer_strings
        unseen_kmer_strings = collections.OrderedDict()
        for kmer in old_unseen_kmer_strings.keys():
            if kmer in initial_kmers:
                unseen_kmer_strings[kmer] = old_unseen_kmer_strings[kmer]

    logger.info('Building subgraphs')
    runner = MccortexRunner(dist=len(unseen_kmer_strings),
                            mem=args.memory,
                            threads=args.cores)
    graph_idx = 0
    while len(unseen_kmer_strings) != 0:
        graph_id = f'g{graph_idx}'
        logger.info(f'Building graph {graph_idx}. Remaining kmers: {len(unseen_kmer_strings)}')

        initial_kmer_string = next(iter(unseen_kmer_strings.keys()))
        subgraph_path = out_dir / f'{graph_id}.ctx'

        runner.build_subgraph(input=input_graph,
                              out=subgraph_path,
                              initial_kmer=initial_kmer_string,
                              graph_id=graph_id)
        remove_subgraph_kmers_from_dict(graph=subgraph_path, dict=unseen_kmer_strings)
        graph_idx += 1
