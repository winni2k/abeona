import attr
import networkx as nx
import collections
from pathlib import Path
from Bio import SeqIO
from cortexpy.graph import traversal
from cortexpy.graph.parser import RandomAccess
from cortexpy.graph.parser.streaming import kmer_string_generator_from_stream
from cortexpy.utils import lexlo
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('abeona.subgraphs')


def load_initial_kmers(input_graph_fh, initial_contigs):
    logger.info('Subsetting graph kmer strings on initial-contig kmers')
    initial_kmers = set()
    first_kmer = next(kmer_string_generator_from_stream(input_graph_fh))
    for rec in SeqIO.parse(initial_contigs, 'fasta'):
        kmer_size = len(first_kmer)
        assert len(rec.seq) >= kmer_size
        for start_idx in range(len(rec.seq) - kmer_size + 1):
            initial_kmers.add(str(lexlo(rec.seq[start_idx:start_idx + kmer_size])))
    logger.info(f'Loaded {len(initial_kmers)} initial kmers')
    ra = RandomAccess(input_graph_fh)
    kmer_strings = KmerStringTracker(
        unseen_generator=sorted(kmer for kmer in initial_kmers if kmer in ra))
    logger.info(
        f'Keeping {len(kmer_strings.unseen_generator)} out of {len(initial_kmers)} initial kmers')
    return kmer_strings


@attr.s(slots=True)
class KmerStringTracker(object):
    unseen_generator = attr.ib(attr.Factory(list))
    seen = attr.ib(attr.Factory(set))

    def add_seen(self, string):
        self.seen.add(string)

    def strings_not_seen(self):
        for string in self.unseen_generator:
            if string not in self.seen:
                yield string


def main(args):
    out_dir = Path(args.out_dir)
    out_dir.mkdir(exist_ok=True)

    input_graph = Path(args.graph)
    if not input_graph.is_file():
        raise Exception(f'Input cortex graph ({input_graph}) does not exist')

    with open(input_graph, 'rb') as first_input_graph_fh:
        if args.initial_contigs:
            kmer_strings = load_initial_kmers(first_input_graph_fh, args.initial_contigs)
        else:
            kmer_strings = KmerStringTracker(
                unseen_generator=kmer_string_generator_from_stream(first_input_graph_fh))
        graph_idx = 0
        with open(input_graph, 'rb') as fh:
            ra_parser = RandomAccess(fh)
            for initial_kmer_string in kmer_strings.strings_not_seen():
                graph_id = f'g{graph_idx}'
                subgraph_path = out_dir / f'{graph_id}.traverse.pickle'
                # logger.info(f'Remaining kmers: {tot_recs}')
                logger.info(f'Building graph {graph_idx} and writing to: {subgraph_path}')

                engine = traversal.Engine(
                    ra_parser,
                    orientation=traversal.constants.EngineTraversalOrientation.both,
                    max_nodes=None,
                    logging_interval=90
                )
                engine.traverse_from(initial_kmer_string)
                logger.info(f'Found subgraph with {len(engine.graph)} kmers')
                with open(subgraph_path, 'wb') as out_fh:
                    nx.write_gpickle(engine.graph, out_fh)
                for node in engine.graph:
                    kmer_strings.add_seen(lexlo(node))
                graph_idx += 1
    logger.info('No kmers remaining')
