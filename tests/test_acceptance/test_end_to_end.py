import gzip
import itertools
from pathlib import Path

import attr
import networkx as nx
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from cortexpy.test.expectation import KmerGraphExpectation

import abeona.__main__
from cortexpy.graph.parser.streaming import kmer_list_generator_from_stream


def get_graph_stream_iterator(file_handle):
    while True:
        try:
            yield nx.read_gpickle(file_handle)
        except EOFError:
            break


@attr.s(slots=True)
class AbeonaExpectation(object):
    out_dir = attr.ib()
    out_graph = attr.ib(init=False)
    out_clean = attr.ib(init=False)

    def __attrs_post_init__(self):
        self.out_dir = Path(self.out_dir)
        self.out_graph = self.out_dir / 'cortex_graph' / 'full.ctx'
        self.out_clean = self.out_graph.with_suffix('.clean.ctx')

    def has_subgraph(self, sg_id):
        return AbeonaSubgraphExpectation(self.out_dir, sg_id)

    def has_n_subgraphs(self, n):
        assert n == len(list((Path(self.out_dir) / 'cortex_subgraphs').glob('g*.ctx')))
        return self

    def has_out_graph_with_kmers(self, *kmers):
        with open(self.out_graph, 'rb') as fh:
            output_kmers = [''.join(k) for k in kmer_list_generator_from_stream(fh)]
        assert set(kmers) == set(output_kmers)
        return self

    def has_out_clean_graph_with_kmers(self, *kmers):
        with open(self.out_clean, 'rb') as fh:
            output_kmers = [''.join(k) for k in kmer_list_generator_from_stream(fh)]
        assert set(kmers) == set(output_kmers)
        return self


@attr.s(slots=True)
class AbeonaSubgraphExpectation(object):
    out_dir = attr.ib()
    sg_id = attr.ib()

    def has_cortex_graph_with_kmers(self, *kmers):
        subgraph = self.out_dir / 'cortex_subgraphs' / f'g{self.sg_id}.ctx'
        assert subgraph.is_file()
        output_kmers = [''.join(k) for k in
                        kmer_list_generator_from_stream(open(subgraph, 'rb'))]
        assert set(kmers) == set(output_kmers)
        return self

    def has_traversal(self):
        subgraph = self.out_dir / 'traversals' / f'g{self.sg_id}.traverse.pickle'
        assert subgraph.is_file()
        graphs = list(get_graph_stream_iterator(open(subgraph, 'rb')))
        assert 1 == len(graphs)
        return KmerGraphExpectation(graphs[0])

    def _has_seqs(self, *seq_strings, dir_name):
        transcripts = self.out_dir / dir_name / f'g{self.sg_id}.transcripts.fa.gz'
        assert transcripts.is_file()
        with gzip.open(str(transcripts), 'rt') as fh:
            seqs = list(SeqIO.parse(fh, 'fasta'))
        assert set(seq_strings) == {str(s.seq) for s in seqs}

    def has_candidate_transcripts(self, *seq_strings):
        self._has_seqs(*seq_strings, dir_name='candidate_transcripts')
        return self

    def has_transcripts(self, *seq_strings):
        self._has_seqs(*seq_strings, dir_name='transcripts')


def traverses_three_subgraphs_into_two_transcripts(tmpdir):
    # given
    kmer_size = 3
    sequences = [('AAAT', 4), ('AAAC', 4)]
    sequences = list(itertools.chain(*[[seq for _ in range(num)] for seq, num in sequences]))
    seq_recs = [SeqRecord(Seq(seq), id=str(idx),
                          letter_annotations={"phred_quality": [40 for _ in range(len(seq))]}) for
                idx, seq in enumerate(sequences)]
    input_fastq = tmpdir / 'single.fq'
    with open(input_fastq, 'w') as fh:
        SeqIO.write(seq_recs, fh, 'fastq')

    expected_sequences = ['AAAT', 'AAAC']
    expected_kmers = [{'AAA', 'AAT', 'AAC'}]
    expected_aligned_kmers = [{'AAA', 'AAT'}, {'ATC', 'TCC'}]

    out_dir = Path(tmpdir) / 'abeona'
    out_graph = out_dir / 'cortex_graph' / 'full.ctx'
    out_clean = out_graph.with_suffix('.clean.ctx')

    # when
    abeona.__main__.main(str(c) for c in ['assemble',
                                          '--fastx-single', input_fastq,
                                          '--kallisto-fastx-single', input_fastq,
                                          '--kallisto-fragment-length', 4,
                                          '--kallisto-sd', 0.01,
                                          '--out-dir', out_dir,
                                          '--kmer-size', kmer_size,
                                          '--min-tip-length', 0,
                                          '--quiet'])

    # then
    with open(out_graph, 'rb') as fh:
        output_kmers = [''.join(k) for k in kmer_list_generator_from_stream(fh)]
    assert {'AAA', 'AAT', 'ATC', 'GGA', 'CCC', 'CCG'} == set(output_kmers)

    with open(out_clean, 'rb') as fh:
        output_kmers = [''.join(k) for k in kmer_list_generator_from_stream(fh)]
    assert {'AAA', 'AAT', 'ATC', 'GGA'} == set(output_kmers)

    # cortex subgraphs
    for sg_id in range(2):
        subgraph = out_dir / 'cortex_subgraphs' / f'g{sg_id}.ctx'
        assert subgraph.is_file()
        output_kmers = [''.join(k) for k in kmer_list_generator_from_stream(open(subgraph, 'rb'))]
        assert expected_kmers[sg_id] == set(output_kmers)

    # traversals
    for sg_id in range(2):
        subgraph = out_dir / 'traversals' / f'g{sg_id}.traverse.pickle'
        assert subgraph.is_file()
        graphs = list(get_graph_stream_iterator(open(subgraph, 'rb')))
        assert 1 == len(graphs)
        expect = KmerGraphExpectation(graphs[0])
        expect.has_nodes(*expected_aligned_kmers[sg_id])

    # transcripts
    for sg_id in range(2):
        transcripts = out_dir / 'transcripts' / f'g{sg_id}.transcripts.fa.gz'
        assert transcripts.is_file()
        with gzip.open(str(transcripts), 'rt') as fh:
            seqs = list(SeqIO.parse(fh, 'fasta'))
        assert 1 == len(seqs)
        assert expected_sequences[sg_id] == str(seqs[0].seq)


def test_traverses_two_subgraphs_into_two_transcripts(tmpdir):
    # given
    kmer_size = 3
    sequences = ['AAAT', 'ATCC']
    seq_recs = [SeqRecord(Seq(seq), id=str(idx),
                          letter_annotations={"phred_quality": [40 for _ in range(len(seq))]}) for
                idx, seq in enumerate(sequences)]
    input_fastq = tmpdir / 'single.fq'
    with open(input_fastq, 'w') as fh:
        SeqIO.write(seq_recs, fh, 'fastq')

    expected_sequences = [['AAAT'], ['ATCC']]
    expected_kmers = [{'AAA', 'AAT'}, {'ATC', 'GGA'}]
    all_expected_kmers = list(itertools.chain(*expected_kmers))
    expected_aligned_kmers = [{'AAA', 'AAT'}, {'ATC', 'TCC'}]

    out_dir = Path(tmpdir) / 'abeona'

    # when
    abeona.__main__.main(str(c) for c in ['assemble',
                                          '--fastx-single', input_fastq,
                                          '--kallisto-fastx-single', input_fastq,
                                          '--kallisto-fragment-length', 3,
                                          '--kallisto-sd', 0.1,
                                          '--bootstrap-samples', 100,
                                          '--out-dir', out_dir,
                                          '--kmer-size', kmer_size,
                                          '--min-tip-length', 0,
                                          '--min-unitig-coverage', 0, ])

    # then
    expect = AbeonaExpectation(out_dir)
    expect.has_out_graph_with_kmers(*all_expected_kmers)
    expect.has_out_clean_graph_with_kmers(*all_expected_kmers)

    for sg_id in range(2):
        sg_expect = expect.has_subgraph(sg_id)
        sg_expect.has_cortex_graph_with_kmers(*expected_kmers[sg_id]) \
            .has_traversal() \
            .has_nodes(*expected_aligned_kmers[sg_id])
        sg_expect.has_candidate_transcripts(*expected_sequences[sg_id])
        sg_expect.has_transcripts(*expected_sequences[sg_id])


def test_when_pruning_traverses_two_subgraphs_into_two_transcripts(tmpdir):
    # given
    kmer_size = 3
    min_tip_length = 2
    sequences = ['AAATG', 'AAAC', 'ATCCC', 'ATCG']
    seq_recs = [SeqRecord(Seq(seq), id=str(idx),
                          letter_annotations={"phred_quality": [40 for _ in range(len(seq))]}) for
                idx, seq in enumerate(sequences)]
    input_fastq = tmpdir / 'single.fq'
    with open(input_fastq, 'w') as fh:
        SeqIO.write(seq_recs, fh, 'fastq')

    expected_sequences = [['AAATG'], ['ATCCC']]
    expected_kmers = [{'AAA', 'AAT', 'AAC', 'ATG'}, {'ATC', 'GGA', 'CCC', 'CGA'}]
    all_expected_kmers = list(itertools.chain(*expected_kmers))
    expected_aligned_kmers = [{'AAA', 'AAT', 'ATG', 'AAC'}, {'ATC', 'TCC', 'CCC', 'TCG'}]

    out_dir = Path(tmpdir) / 'abeona'

    # when
    abeona.__main__.main(str(c) for c in ['assemble',
                                          '--fastx-single', input_fastq,
                                          '--kallisto-fastx-single', input_fastq,
                                          '--kallisto-fragment-length', 3,
                                          '--kallisto-sd', 0.1,
                                          '--bootstrap-samples', 100,
                                          '--out-dir', out_dir,
                                          '--kmer-size', kmer_size,
                                          '--min-tip-length', min_tip_length,
                                          '--min-unitig-coverage', 0,
                                          '--quiet'])

    # then
    expect = AbeonaExpectation(out_dir)
    expect.has_out_graph_with_kmers(*all_expected_kmers)
    expect.has_out_clean_graph_with_kmers(*all_expected_kmers)

    for sg_id in range(2):
        sg_expect = expect.has_subgraph(sg_id)
        sg_expect.has_cortex_graph_with_kmers(*expected_kmers[sg_id]) \
            .has_traversal() \
            .has_nodes(*expected_aligned_kmers[sg_id])
        sg_expect.has_candidate_transcripts(*expected_sequences[sg_id])
        sg_expect.has_transcripts(*expected_sequences[sg_id])


def test_when_filtering_traverses_three_subgraphs_into_two_transcripts(tmpdir):
    # given
    kmer_size = 3
    sequences = [('AAAT', 4), ('ATCC', 4), ('CCCG', 3)]
    sequences = list(itertools.chain(*[[seq for _ in range(num)] for seq, num in sequences]))
    seq_recs = [SeqRecord(Seq(seq), id=str(idx),
                          letter_annotations={"phred_quality": [40 for _ in range(len(seq))]}) for
                idx, seq in enumerate(sequences)]
    input_fastq = tmpdir / 'single.fq'
    with open(input_fastq, 'w') as fh:
        SeqIO.write(seq_recs, fh, 'fastq')

    expected_sequences = [['AAAT'], ['ATCC']]

    all_expected_kmers = {'AAA', 'AAT', 'ATC', 'GGA', 'CCC', 'CCG'}
    expected_clean_kmers = [{'AAA', 'AAT'}, {'ATC', 'GGA'}]
    all_expected_clean_kmers = list(itertools.chain(*expected_clean_kmers))
    expected_aligned_kmers = [{'AAA', 'AAT'}, {'ATC', 'TCC'}]

    out_dir = Path(tmpdir) / 'abeona'

    # when
    abeona.__main__.main(str(c) for c in ['assemble',
                                          '--fastx-single', input_fastq,
                                          '--kallisto-fastx-single', input_fastq,
                                          '--kallisto-fragment-length', 3,
                                          '--kallisto-sd', 0.1,
                                          '--bootstrap-samples', 100,
                                          '--out-dir', out_dir,
                                          '--kmer-size', kmer_size,
                                          '--min-tip-length', 0])

    # then
    expect = AbeonaExpectation(out_dir)
    expect.has_out_graph_with_kmers(*all_expected_kmers)
    expect.has_out_clean_graph_with_kmers(*all_expected_clean_kmers)

    for sg_id in range(2):
        sg_expect = expect.has_subgraph(sg_id)
        sg_expect.has_cortex_graph_with_kmers(*expected_clean_kmers[sg_id]) \
            .has_traversal() \
            .has_nodes(*expected_aligned_kmers[sg_id])
        sg_expect.has_candidate_transcripts(*expected_sequences[sg_id])
        sg_expect.has_transcripts(*expected_sequences[sg_id])


class TestInitialSeqsFasta(object):
    def test_traverses_two_subgraphs_into_single_transcript(self, tmpdir):
        # given
        kmer_size = 3
        sequences = ['AAAT', 'ATCC']
        initial_seqs = ['AAAT']

        inital_contigs_fasta = tmpdir / 'initial-contigs.fa'
        SeqIO.write([SeqRecord(Seq(seq), id=str(idx)) for idx, seq in enumerate(initial_seqs)],
                    str(inital_contigs_fasta),
                    'fasta')

        input_fastq = tmpdir / 'single.fq'
        seq_recs = [SeqRecord(Seq(seq), id=str(idx),
                              letter_annotations={"phred_quality": [40 for _ in range(len(seq))]})
                    for
                    idx, seq in enumerate(sequences)]
        with open(input_fastq, 'w') as fh:
            SeqIO.write(seq_recs, fh, 'fastq')

        expected_kmers = [{'AAA', 'AAT'}, {'ATC', 'GGA'}]
        all_expected_kmers = list(itertools.chain(*expected_kmers))

        expected_subgraph_sequences = [['AAAT']]
        expected_subgraph_kmers = [{'AAA', 'AAT'}]
        expected_subgraph_aligned_kmers = [{'AAA', 'AAT'}, {'ATC', 'TCC'}]

        out_dir = Path(tmpdir) / 'abeona'

        # when
        abeona.__main__.main(str(c) for c in ['assemble',
                                              '--initial-contigs', inital_contigs_fasta,
                                              '--fastx-single', input_fastq,
                                              '--kallisto-fastx-single', input_fastq,
                                              '--kallisto-fragment-length', 3,
                                              '--kallisto-sd', 0.1,
                                              '--bootstrap-samples', 100,
                                              '--out-dir', out_dir,
                                              '--kmer-size', kmer_size,
                                              '--min-tip-length', 0,
                                              '--min-unitig-coverage', 0, ])

        # then
        expect = AbeonaExpectation(out_dir)
        expect.has_out_graph_with_kmers(*all_expected_kmers)
        expect.has_out_clean_graph_with_kmers(*all_expected_kmers)

        for sg_id in range(1):
            sg_expect = expect.has_subgraph(sg_id)
            sg_expect.has_cortex_graph_with_kmers(*expected_subgraph_kmers[sg_id]) \
                .has_traversal() \
                .has_nodes(*expected_subgraph_aligned_kmers[sg_id])
            sg_expect.has_candidate_transcripts(*expected_subgraph_sequences[sg_id])
            sg_expect.has_transcripts(*expected_subgraph_sequences[sg_id])
        expect.has_n_subgraphs(1)
