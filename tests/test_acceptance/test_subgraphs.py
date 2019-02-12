import collections

from hypothesis import given, strategies as strat

from abeona.test.drivers import SubgraphTestDriver

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
    d = SubgraphTestDriver(tmpdir)
    for seq in seqs:
        d.with_dna_sequence(seq.seq)

    # when
    expect = d.run()

    # then
    for seq in seqs:
        expect.has_graph_with_kmers(*seq.kmers).has_meta_info('n_junctions', "0")
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


def test_counts_one_junction_of_subgraphs(tmpdir):
    # given
    d = SubgraphTestDriver(tmpdir)
    d.with_dna_sequence('AAAT')
    d.with_dna_sequence('AAAC')

    # when
    expect = d.run()

    # then
    expect.has_graph_with_kmers('AAA', 'AAT', 'AAC').has_meta_info('n_junctions', "1")
    expect.has_n_graphs(1)


def test_counts_two_junctions_of_subgraphs(tmpdir):
    # given
    d = SubgraphTestDriver(tmpdir)
    d.with_dna_sequence('AAAT')
    d.with_dna_sequence('AAAC')
    d.with_dna_sequence('AACA')
    d.with_dna_sequence('AACT')

    # when
    expect = d.run()

    # then
    expect \
        .has_graph_with_kmers('AAA', 'AAT', 'AAC', 'ACA', 'ACT') \
        .has_meta_info('n_junctions', "2")
    expect.has_n_graphs(1)
