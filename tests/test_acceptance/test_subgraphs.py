import collections
from cortexpy.test import builder
from hypothesis import given, strategies as strat

import abeona.__main__
from tests.expectation.mccortex import Graphs, Graph

Seq = collections.namedtuple('Seq', ['seq', 'order', 'kmers'])


@given(strat.sets(
    strat.sampled_from([
        Seq('AAAT', 0, ('AAA', 'AAT')),
        Seq('ATCC', 1, ('ATC', 'GGA')),
        Seq('CCCG', 2, ('CCC', 'CCG')),
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
    graphs = [out_dir / f'g{i}.ctx' for i in range(len(seqs))]

    # when
    abeona.__main__.main(['subgraphs',str(graph), str(out_dir)])
    expect = Graphs(graphs)

    # then
    for seq in seqs:
        expect.has_graph_with_kmers(*seq.kmers)
    expect.has_n_graphs(len(seqs))

    for graph, seq in zip(graphs, sorted(seqs, key=lambda s: s.order)):
        Graph(graph).has_kmer_strings(*seq.kmers)
