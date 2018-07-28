import attr
import networkx as nx
from pathlib import Path
from Bio import SeqIO
from cortexpy.graph import traversal
from cortexpy.graph.parser import RandomAccess
from cortexpy.graph.serializer.kmer import dump_colored_de_bruijn_graph_to_cortex
from cortexpy.utils import lexlo
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('abeona.subgraphs')


def load_initial_kmers(random_access, initial_contigs):
    logger.info('Subsetting graph kmer strings on initial-contig kmers')
    initial_kmers = set()
    first_kmer = next(iter(random_access))
    for rec in SeqIO.parse(initial_contigs, 'fasta'):
        kmer_size = len(first_kmer)
        assert len(rec.seq) >= kmer_size
        for start_idx in range(len(rec.seq) - kmer_size + 1):
            initial_kmers.add(str(lexlo(rec.seq[start_idx:start_idx + kmer_size])))
    logger.info(f'Loaded {len(initial_kmers)} initial kmers')
    k_string_tracker = KmerStringTracker(
        unseen=sorted(kmer for kmer in initial_kmers if kmer in random_access))
    logger.info(
        f'Keeping {len(k_string_tracker.unseen)} out of {len(initial_kmers)} initial kmers')
    return k_string_tracker


@attr.s(slots=True)
class KmerStringTracker(object):
    unseen = attr.ib(attr.Factory(list))
    seen = attr.ib(attr.Factory(set))

    def add_seen(self, string):
        self.seen.add(string)

    def strings_not_seen(self):
        for string in self.unseen:
            if string not in self.seen:
                yield string


def main(args):
    out_dir = Path(args.out_dir)
    out_dir.mkdir(exist_ok=True)

    input_graph = Path(args.graph)
    if not input_graph.is_file():
        raise Exception(f'Input cortex graph ({input_graph}) does not exist')

    with open(input_graph, 'rb') as first_input_graph_fh:
        ra = RandomAccess(first_input_graph_fh, kmer_cache_size=0)
        if args.initial_contigs:
            kstring_tracker = load_initial_kmers(ra, args.initial_contigs)
        else:
            kstring_tracker = KmerStringTracker(unseen=ra)
        graph_idx = 0
        with open(input_graph, 'rb') as fh:
            ra_parser = RandomAccess(fh, kmer_cache_size=0)
            for initial_kmer_string in kstring_tracker.strings_not_seen():
                graph_id = f'g{graph_idx}'
                subgraph_path = out_dir / f'{graph_id}.traverse.ctx'
                logger.info(f'Building graph {graph_idx} and writing to: {subgraph_path}')

                engine = traversal.Engine(
                    ra_parser,
                    orientation=traversal.constants.EngineTraversalOrientation.both,
                    max_nodes=None,
                    logging_interval=90
                )
                engine.traverse_from(str(initial_kmer_string))
                with open(subgraph_path, 'wb') as out_fh:
                    dump_colored_de_bruijn_graph_to_cortex(engine.graph, out_fh)
                for node in engine.graph:
                    kstring_tracker.add_seen(lexlo(node))
                logger.info(
                    f'Found subgraph with {len(engine.graph)} kmers - at most {len(ra_parser) - len(kstring_tracker.seen)} kmers left')
                graph_idx += 1
    logger.info('No kmers remaining')
