import collections

from hypothesis import given, strategies as strat

from abeona.test.drivers import ReadsTestDriver

SeqTup = collections.namedtuple('SeqTup', ['seq', 'kmers'])


# @given(strat.sets(
#     strat.sampled_from([
#         SeqTup('AAAT', ('AAA', 'AAT')),
#         SeqTup('ATCC', ('ATC', 'GGA')),
#         SeqTup('CCCG', ('CCC', 'CCG')),
#     ]),
#     min_size=1,
#     max_size=3))
# def test_decompose_and_assign_reads_to_graphs(tmpdir, seqs):
#     # given
#     d = ReadsTestDriver(tmpdir)
#     for seq in seqs:
#         d.with_dna_sequence(seq.seq)
#
#     # when
#     expect = d.run()
#
#     # then
#     for seq in seqs:
#         expect.has_graph_with_kmers(*seq.kmers)
#         expect.has_
#     expect.has_n_graphs(len(seqs))
#


