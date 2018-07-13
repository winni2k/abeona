import collections
import shutil
from pathlib import Path

import attr
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from cortexpy.test import builder
from hypothesis import given, strategies as strat

import abeona.__main__
from tests.expectation.mccortex import Traversals

SeqTup = collections.namedtuple('SeqTup', ['seq', 'kmers'])


@attr.s(slots=True)
class SubgraphTestDriver(object):
    tmpdir = attr.ib()
    builder = attr.ib(attr.Factory(builder.Mccortex))
    initial_contigs = attr.ib(None)

    def __getattr__(self, item):
        if item in ['with_kmer_size', 'with_dna_sequence']:
            return getattr(self.builder, item)
        else:
            raise ValueError(f'Could not find {item}')

    def with_initial_contigs(self, *contigs):
        self.initial_contigs = contigs
        return self

    def run(self):
        shutil.rmtree(self.tmpdir)
        self.tmpdir.mkdir()
        graph = self.builder.build(self.tmpdir)
        out_dir = self.tmpdir / 'abeona_output'
        command = ['subgraphs', graph, out_dir]

        # graphs = [out_dir / f'g{i}.traverse.ctx' for i in range(self.builder.num_sequences)]

        if self.initial_contigs:
            initial_contigs = self.tmpdir / 'initial-kmers.fa'
            records = [
                SeqRecord(Seq(contig), id=f'initial_contig_{idx}') for idx, contig in
                enumerate(self.initial_contigs)
            ]
            SeqIO.write(records, str(initial_contigs), 'fasta')
            command += ['--initial-contigs', initial_contigs]

        abeona.__main__.main([str(arg) for arg in command])

        graphs = list(Path(out_dir).glob('g*.traverse.ctx'))

        return Traversals(graphs)


@given(strat.sets(
    strat.sampled_from([
        SeqTup('AAAT', ('AAA', 'AAT')),
        SeqTup('ATCC', ('ATC', 'GGA')),
        SeqTup('CCCG', ('CCC', 'CCG')),
    ]),
    min_size=1,
    max_size=3))
def test_decompose(tmpdir, seqs):
    # given
    d = SubgraphTestDriver(tmpdir)
    for seq in seqs:
        d.with_dna_sequence(seq.seq)

    # when
    expect = d.run()

    # then
    for seq in seqs:
        expect.has_graph_with_kmers(*seq.kmers)
    expect.has_n_graphs(len(seqs))


def test_allows_initial_kmer_not_in_seqs(tmpdir):
    # given
    d = SubgraphTestDriver(tmpdir)
    seqs = [
        SeqTup('AAAT', ('AAA', 'AAT')),
        SeqTup('ATCC', ('ATC', 'GGA')),
        SeqTup('CCCG', ('CCC', 'CCG')),
    ]
    for seq in seqs:
        d.with_dna_sequence(seq.seq)
    d.with_initial_contigs('AAAC')

    # when
    expect = d.run()

    # then
    expect.has_graph_with_kmers(*seqs[0].kmers)
    expect.has_n_graphs(1)
