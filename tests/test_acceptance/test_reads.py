import attr
from cortexpy.utils import lexlo
from hypothesis import given, strategies as strat

from abeona.test.drivers import ReadsTestDriver


@attr.s(slots=True)
class SeqTup:
    seq = attr.ib()
    seq2 = attr.ib(None)
    kmer_size = attr.ib(3)

    def kmers(self):
        kmers = set()
        seqs = [self.seq]
        if self.seq2 is not None:
            seqs.append(self.seq2)
        for seq in seqs:
            for start in range(len(seq) - self.kmer_size + 1):
                kmers.add(seq[start:start + self.kmer_size])

        return kmers

    def lexlo_kmers(self):
        return {lexlo(k) for k in self.kmers()}


@given(strat.sets(
    strat.sampled_from(['AAAT', 'ATCC', 'CCCG']),
    min_size=1,
    max_size=3))
def test_with_single_end_reads_decompose_and_assign_reads_to_graphs(tmpdir, seqs):
    # given
    seqs = [SeqTup(s) for s in seqs]
    d = ReadsTestDriver(tmpdir)
    for seq in seqs:
        d.with_dna_sequence(seq.seq)

    # when
    expect = d.run()

    # then
    for seq in seqs:
        expect.has_graph_with_kmers_and_mapped_reads(seq.lexlo_kmers(), [seq.seq])
    expect.has_n_graphs(len(seqs))


@given(strat.sets(
    strat.sampled_from([
        ('AAAT', 'AATA'),
        ('ATCC', 'TCCA'),
        ('CCCG', 'CCGA'),
    ]),
    min_size=1,
    max_size=3))
def test_with_paired_end_reads_decompose_and_assign_reads_to_graphs(tmpdir, seqs):
    # given
    seqs = [SeqTup(seq=s[0], seq2=s[1]) for s in seqs]
    d = ReadsTestDriver(tmpdir)
    for seq in seqs:
        d.with_dna_sequence_pair(seq.seq, seq.seq2)

    # when
    expect = d.run()

    # then
    for seq in seqs:
        expect.has_graph_with_kmers_and_mapped_reads(seq.lexlo_kmers(), [seq.seq], [seq.seq2])
    expect.has_n_graphs(len(seqs))
