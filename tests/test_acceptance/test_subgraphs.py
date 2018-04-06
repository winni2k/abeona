import collections

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from cortexpy.test import builder
from hypothesis import given, strategies as strat

import abeona.__main__
from tests.expectation.mccortex import Traversals

SeqTup = collections.namedtuple('SeqTup', ['seq', 'kmers'])


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
    b = builder.Mccortex(mccortex_bin='mccortex')
    for seq in seqs:
        b.with_dna_sequence(seq.seq)
    graph = b.build(tmpdir)

    out_dir = tmpdir / 'abeona_output'
    graphs = [out_dir / f'g{i}.traverse.pickle' for i in range(len(seqs))]

    # when
    abeona.__main__.main(['subgraphs', str(graph), str(out_dir)])
    expect = Traversals(graphs)

    # then
    for seq in seqs:
        expect.has_graph_with_kmers(*seq.kmers)
    expect.has_n_graphs(len(seqs))


def test_allows_initial_kmer_not_in_seqs(tmpdir):
    # given
    b = builder.Mccortex(mccortex_bin='mccortex')
    seqs = [
        SeqTup('AAAT', ('AAA', 'AAT')),
        SeqTup('ATCC', ('ATC', 'GGA')),
        SeqTup('CCCG', ('CCC', 'CCG')),
    ]
    for seq in seqs:
        b.with_dna_sequence(seq.seq)
    graph = b.build(tmpdir)

    initial_contigs = tmpdir / 'initial-kmers.fa'
    SeqIO.write([SeqRecord(Seq('AAAC'), id='initial_contig')], str(initial_contigs), 'fasta')

    out_dir = tmpdir / 'abeona_output'
    graphs = [out_dir / f'g{i}.traverse.pickle' for i in range(len(seqs))]

    # when
    abeona.__main__.main(f'subgraphs {graph} {out_dir} --initial-contigs {initial_contigs}'.split())
    expect = Traversals([graphs[0]])

    # then
    expect.has_graph_with_kmers(*seqs[0].kmers)
    expect.has_n_graphs(1)
