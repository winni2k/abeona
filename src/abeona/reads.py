import gzip
from logging import getLogger
from pathlib import Path

import attr
from Bio import SeqIO
from cortexpy.graph.parser.random_access import RandomAccess

from .utils import get_maybe_gzipped_file_handle

logger = getLogger('abeona.reads')


@attr.s(slots=True)
class GraphData:
    prefix = attr.ib()
    graph_path = attr.ib()
    format = attr.ib()
    is_paired = attr.ib()
    buffer1 = attr.ib(attr.Factory(list))
    buffer2 = attr.ib(attr.Factory(list))
    _fh = attr.ib(None)

    @property
    def fh(self):
        if self._fh is None:
            Path(self.prefix).parent.mkdir(exist_ok=True)
            self._fh = [gzip.open(self.prefix + f'.1.{self.format}.gz', 'wt')]
            if self.is_paired:
                self._fh.append(gzip.open(self.prefix + f'.2.{self.format}.gz', 'wt'))
        return self._fh

    def store_record_pair(self, *records):
        self.buffer1.append(records[0])
        if self.is_paired:
            self.buffer2.append(records[1])

    def write_record_pair(self, *records):
        """Accepts one or two reads. Figures out if one or two reads are expected."""
        self.fh[0].write(records[0].format(self.format))
        if self.is_paired:
            self.fh[1].write(records[1].format(self.format))

    def dump(self):
        SeqIO.write(self.buffer1, self.fh[0], self.format)
        self.fh[0].close()
        if self.is_paired:
            SeqIO.write(self.buffer2, self.fh[1], self.format)
            self.fh[1].close()
        self._fh = None


def main(args):
    is_paired = args.reverse is not None
    graphs = []
    with open(args.graph_list) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            prefix, graph = line.rstrip().split('\t')
            graphs.append(GraphData(prefix=prefix,
                                    graph_path=graph,
                                    format=args.format,
                                    is_paired=is_paired))
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

    assert kmer_size, "Could not find any graphs"

    reads = [SeqIO.parse(get_maybe_gzipped_file_handle(args.forward, 'rt'), args.format)]
    if is_paired:
        reads.append(SeqIO.parse(get_maybe_gzipped_file_handle(args.reverse, 'rt'), args.format))
    for recs in zip(*reads):
        for start in range(len(recs[0]) - kmer_size + 1):
            kmer = recs[0].seq[start:(start + kmer_size)]
            if kmer in kmers:
                graphs[kmers[kmer]].store_record_pair(*recs)
                break

    for graph in graphs:
        graph.dump()
