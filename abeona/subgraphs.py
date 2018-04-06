import collections

import attr
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from cortexpy.graph import traversal
from cortexpy.graph.parser import RandomAccess
from cortexpy.graph.parser.streaming import kmer_string_generator_from_stream
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


def remove_subgraph_kmers_from_dict(graph, kmer_dict):
    """Reads kmer strings from graph and removes the strings from dict in place"""
    with open(graph, 'rb') as fh:
        for kmer in kmer_string_generator_from_stream(fh):
            if kmer in kmer_dict:
                del kmer_dict[kmer]


def main(argv):
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

    from cortexpy.utils import lexlo
    import networkx as nx

    out_dir = Path(args.out_dir)
    out_dir.mkdir(exist_ok=True)

    input_graph = Path(args.graph)
    if not input_graph.is_file():
        raise Exception(f'Input cortex graph ({input_graph}) does not exist')

    if args.initial_contigs:
        logger.info('Subsetting graph kmer strings on initial-contig kmers')
        initial_kmers = set()
        with open(input_graph, 'rb') as fh:
            first_kmer = next(kmer_string_generator_from_stream(fh))
            for rec in SeqIO.parse(args.initial_contigs, 'fasta'):
                kmer_size = len(first_kmer)
                assert len(rec.seq) >= kmer_size
                for start_idx in range(len(rec.seq) - kmer_size + 1):
                    initial_kmers.add(str(lexlo(rec.seq[start_idx:start_idx + kmer_size])))
            logger.info(f'Loaded {len(initial_kmers)} initial kmers')
            ra = RandomAccess(fh)
            unseen_kmer_strings = collections.OrderedDict.fromkeys(
                sorted(kmer for kmer in initial_kmers if kmer in ra))
        logger.info(f'Keeping {len(unseen_kmer_strings)} out of {len(initial_kmers)} initial kmers')
    else:
        logger.info('Loading graph kmer strings')
        with open(input_graph, 'rb') as fh:
            unseen_kmer_strings = collections.OrderedDict.fromkeys(
                kmer_string_generator_from_stream(fh))
        logger.info(f'Completed loading {len(unseen_kmer_strings)} kmers')

    logger.info(f'Building subgraphs from {len(unseen_kmer_strings)} initial kmers')
    graph_idx = 0
    with open(input_graph, 'rb') as fh:
        ra_parser = RandomAccess(fh)
        while len(unseen_kmer_strings) != 0:
            graph_id = f'g{graph_idx}'
            subgraph_path = out_dir / f'{graph_id}.traverse.pickle'
            logger.info(f'Remaining kmers: {len(unseen_kmer_strings)}')
            logger.info(f'Building graph {graph_idx} and writing to: {subgraph_path}')

            initial_kmer_string = next(iter(unseen_kmer_strings.keys()))
            engine = traversal.Engine(
                ra_parser,
                orientation=traversal.constants.EngineTraversalOrientation.both,
                max_nodes=None,
                logging_interval=90
            )
            engine.traverse_from(initial_kmer_string)
            with open(subgraph_path, 'wb') as out_fh:
                nx.write_gpickle(engine.graph, out_fh)
            for node in engine.graph:
                node = lexlo(node)
                if node in unseen_kmer_strings.keys():
                    del unseen_kmer_strings[node]
            graph_idx += 1
    logger.info('No kmers remaining')
