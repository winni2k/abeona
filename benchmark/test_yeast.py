import itertools
import pytest

import abeona.subgraphs

CHROM_GRAPH_230KB = 'fixtures/yeast/NC_001133.9.ctx'
CHROM_GRAPH_1KB = 'fixtures/yeast/NC_001133.9.1kbp.ctx'
CHROM_GRAPH_4KB = 'fixtures/yeast/NC_001133.9.4kbp.ctx'


def run_abeona_subgraphs(graph, out_dir):
    cmd = ['-m', '3', '--cores', '2', str(graph), str(out_dir)]
    abeona.subgraphs.main(cmd)


GRAPHS = {
    '230k': CHROM_GRAPH_230KB,
    '1kb': CHROM_GRAPH_1KB,
    '4kb': CHROM_GRAPH_4KB
}

FUNCS = {
    'subgraphs': run_abeona_subgraphs
}


@pytest.mark.parametrize('graph_size,func_type',
                         itertools.product(GRAPHS.keys(), FUNCS.keys()))
def test_yeast(benchmark, tmpdir, graph_size, func_type):
    graph = GRAPHS[graph_size]
    func = FUNCS[func_type]
    benchmark(func, graph, tmpdir)
