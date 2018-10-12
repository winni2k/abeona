import gzip
from logging import getLogger
from pathlib import Path

import attr
from Bio import SeqIO
from cortexpy.graph.parser.random_access import RandomAccess

logger = getLogger('abeona.reads')


@attr.s(slots=True)
class GraphData:
    prefix = attr.ib()
    graph_path = attr.ib()
    format = attr.ib()
    _fh = attr.ib(None)

    @property
    def fh(self):
        if self._fh is None:
            Path(self.prefix).parent.mkdir(exist_ok=True)
            self._fh = gzip.open(self.prefix + f'.{self.format}.gz', 'wt')
        return self._fh

    def write_record(self, record):
        self.fh.write(record.format(self.format))


def main(args):
    graphs = []
    with open(args.graph_list) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            prefix, graph = line.rstrip().split('\t')
            graphs.append(GraphData(prefix=prefix,
                                    graph_path=graph,
                                    format=args.format))
    logger.info(f'Assigning reads to {len(graphs)} graphs')

    kmers = {}
    kmer_size = None
    for idx, graph_data in enumerate(graphs):
        with open(graph_data.graph_path, 'rb') as fh:
            ra = RandomAccess(fh)
            if kmer_size is None:
                kmer_size = ra.kmer_size
                logger.info(f'Detected kmer size: {kmer_size}')
            assert kmer_size == ra.kmer_size
            for kmer in ra:
                kmers[kmer] = idx

    for rec in SeqIO.parse(args.forward, args.format):
        for start in range(len(rec) - kmer_size + 1):
            kmer = rec.seq[start:(start + kmer_size)]
            if kmer in kmers:
                graphs[kmers[kmer]].write_record(rec)
                break
