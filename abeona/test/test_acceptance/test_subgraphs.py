import os
from cortexpy.test import builder

from abeona.test.expectation.mccortex import Graphs


def test_decompose(tmpdir):
    # given
    b = builder.Mccortex(mccortex_bin='mccortex')
    b.with_dna_sequence('AAAT')
    b.with_dna_sequence('CCCG')
    graph = b.build(tmpdir)

    out_dir = tmpdir / 'abeona_output'
    graph1 = out_dir / 'g0.ctx'
    graph2 = out_dir / 'g1.ctx'

    # when
    os.system(f'abeona subgraphs {graph} {out_dir}')
    expect = Graphs([graph1, graph2])

    # then
    expect.has_graph_with_kmers('AAA', 'AAT')
    expect.has_graph_with_kmers('CCC', 'CCG')
    expect.has_n_graphs(2)
