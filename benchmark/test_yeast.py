import pytest
from abeona.cli import main

CHROM_GRAPH_230KB = 'fixtures/yeast/NC_001133.9.ctx'

GRAPHS = {
    str(kbp) + 'kb': 'fixtures/yeast/NC_001133.9.{}kbp.ctx'.format(kbp)
    for kbp in [1, 2, 4, 8]
}


@pytest.mark.parametrize('graph_size', GRAPHS.keys())
def test_subgraphs_of_yeast(benchmark, tmpdir, graph_size):
    graph = GRAPHS[graph_size]
    args = ['abeona', 'subgraphs', graph, str(tmpdir)]
    benchmark.pedantic(main, args=(args,), iterations=2,
                       warmup_rounds=1)
