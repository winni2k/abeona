import gzip
from logging import getLogger
from pathlib import Path

import attr
from cortexpy.graph.parser.random_access import RandomAccess

from .utils import get_maybe_gzipped_file_handle, get_fastx_title_seq_generator

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
    _file1 = attr.ib(init=False)
    _file2 = attr.ib(init=False)

    def __attrs_post_init__(self):
        Path(self.prefix).parent.mkdir(exist_ok=True)
        self._file1 = self.prefix + f'.1.fa.gz'
        path = Path(self._file1)
        if path.is_file():
            path.unlink()
        if self.is_paired:
            self._file2 = self.prefix + f'.2.fa.gz'
            path = Path(self._file2)
            if path.is_file():
                path.unlink()

    @property
    def fh(self):
        if self._fh is None:
            self._fh = [gzip.open(self._file1, 'at')]
            if self.is_paired:
                self._fh.append(gzip.open(self._file2, 'at'))
        return self._fh

    def store_record_pair(self, *records):
        self.buffer1.append(records[0])
        if self.is_paired:
            self.buffer2.append(records[1])

    def write_record_pair(self, *records):
        """Accepts one or two reads. Figures out if one or two reads are expected."""
        self.fh[0].write(''.join(['>', records[0][0], '\n', records[0][1], '\n']))
        if self.is_paired:
            self.fh[1].write(''.join(['>', records[1][0], '\n', records[1][1], '\n']))

    def flush(self):
        buffers = [self.buffer1]
        if self.is_paired:
            buffers.append(self.buffer2)
        for recs in zip(*buffers):
            self.write_record_pair(*recs)
        self.fh[0].close()
        self.buffer1.clear()
        if self.is_paired:
            self.fh[1].close()
            self.buffer2.clear()
        self._fh = None


def main(args):
    record_buffer_size = args.record_buffer_size

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
                assert kmer not in kmers, 'Input subgraphs must be disjoint'
                kmers[kmer] = idx

    assert kmer_size, "Could not find any graphs"

    reads = [get_fastx_title_seq_generator(get_maybe_gzipped_file_handle(args.forward, 'rt'))]
    if is_paired:
        reads.append(
            get_fastx_title_seq_generator(get_maybe_gzipped_file_handle(args.reverse, 'rt')))
    n_records_stored = 0
    for recs in zip(*reads):
        seen_graphs = set()
        for read in (s for t, s in recs):
            for start in range(len(read) - kmer_size + 1):
                kmer = read[start:(start + kmer_size)]
                if kmer in kmers:
                    if record_buffer_size is not None and record_buffer_size >= n_records_stored:
                        for graph in graphs:
                            graph.flush()
                        n_records_stored = 0
                    idx = kmers[kmer]
                    if idx not in seen_graphs:
                        graphs[idx].store_record_pair(*recs)
                        n_records_stored += 1
                        seen_graphs.add(idx)
    for graph in graphs:
        graph.flush()
