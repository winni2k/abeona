import gzip
import itertools
from pathlib import Path

import attr
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from cortexpy.test.expectation import KmerGraphExpectation
from cortexpy.utils import lexlo

import abeona.__main__
from cortexpy.graph.parser.streaming import kmer_string_generator_from_stream, load_cortex_graph


@attr.s(slots=True)
class AbeonaExpectation(object):
    out_dir = attr.ib()
    prune_length = attr.ib(0)
    out_graph = attr.ib(init=False)
    out_clean = attr.ib(init=False)
    subgraphs = attr.ib(init=False)
    traversal_dir = attr.ib(init=False)
    out_clean_tip_pruned = attr.ib(init=False)

    def __attrs_post_init__(self):
        self.out_dir = Path(self.out_dir)
        self.out_graph = self.out_dir / 'cortex_graph' / 'full.ctx'
        self.out_clean = self.out_graph.with_suffix('.clean.ctx')
        self.out_clean_tip_pruned = self.out_clean.with_suffix(
            f'.min_tip_length_{self.prune_length}.ctx')
        self.traversal_dir = self.out_dir / 'traversals'
        self.subgraphs = [self.traversal_dir / f'g{i}.traverse.ctx'
                          for i, _ in
                          enumerate(self.traversal_dir.glob(f'g*.traverse.ctx'))]

    def has_subgraph(self, sg_id):
        return AbeonaSubgraphExpectation(self.out_dir, sg_id)

    def has_subgraph_with_kmers(self, *kmers):
        kmers = set(kmers)
        assert len(kmers) > 0
        for sg_id, sg in enumerate(self.subgraphs):
            with open(sg, 'rb') as fh:
                subgraph_kmers = set(lexlo(n) for n in load_cortex_graph(fh))
            if len(kmers & subgraph_kmers) > 0:
                assert kmers == subgraph_kmers
                return AbeonaSubgraphExpectation(self.out_dir, sg_id)
        assert kmers == set()
        return self

    def has_n_subgraphs(self, n):
        assert n == len(self.subgraphs)
        return self

    def has_out_graph_with_kmers(self, *kmers):
        with open(self.out_graph, 'rb') as fh:
            output_kmers = list(kmer_string_generator_from_stream(fh))
        assert set(kmers) == set(output_kmers)
        return self

    def has_out_clean_graph_with_kmers(self, *kmers):
        with open(self.out_clean, 'rb') as fh:
            output_kmers = list(kmer_string_generator_from_stream(fh))
        assert set(kmers) == set(output_kmers)
        return self


@attr.s(slots=True)
class AbeonaSubgraphExpectation(object):
    out_dir = attr.ib()
    sg_id = attr.ib()

    def has_cortex_graph_with_kmers(self, *kmers):
        subgraph = self.out_dir / 'cortex_subgraphs' / f'g{self.sg_id}.ctx'
        assert subgraph.is_file()
        output_kmers = list(kmer_string_generator_from_stream(open(subgraph, 'rb')))
        assert set(kmers) == set(output_kmers)
        return self

    def has_traversal(self):
        subgraph = self.out_dir / 'traversals' / f'g{self.sg_id}.traverse.ctx'
        assert subgraph.is_file()
        graph = load_cortex_graph(open(subgraph, 'rb'))
        return KmerGraphExpectation(graph)

    def _has_seqs(self, *expected_seq_strings, dir_name):
        transcripts = self.out_dir / dir_name / f'g{self.sg_id}.transcripts.fa.gz'
        assert transcripts.is_file()
        with gzip.open(str(transcripts), 'rt') as fh:
            seqs = list(SeqIO.parse(fh, 'fasta'))
        seq_strings = {str(s.seq) for s in seqs} | {str(s.seq.reverse_complement()) for s in seqs}
        for expected_seq_string in expected_seq_strings:
            assert expected_seq_string in seq_strings
        assert len(set(expected_seq_strings)) == len(seqs)

    def has_candidate_transcripts(self, *seq_strings):
        self._has_seqs(*seq_strings, dir_name='candidate_transcripts')
        return self

    def has_transcripts(self, *seq_strings):
        self._has_seqs(*seq_strings, dir_name='transcripts')


class TestAssemble(object):
    def test_traverses_two_subgraphs_into_two_transcripts(self, tmpdir):
        # given
        kmer_size = 3
        sequences = ['AAAT', 'ATCC']
        seq_recs = [SeqRecord(Seq(seq), id=str(idx),
                              letter_annotations={"phred_quality": [40 for _ in range(len(seq))]})
                    for
                    idx, seq in enumerate(sequences)]
        input_fastq = tmpdir / 'single.fq'
        with open(input_fastq, 'w') as fh:
            SeqIO.write(seq_recs, fh, 'fastq')

        expected_sequences = [['AAAT'], ['ATCC']]
        expected_kmers = [{'AAA', 'AAT'}, {'ATC', 'GGA'}]
        all_expected_kmers = list(itertools.chain(*expected_kmers))

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
            sg_expect = expect.has_subgraph_with_kmers(*expected_kmers[sg_id])
            sg_expect.has_traversal().has_nodes(*expected_kmers[sg_id])
            sg_expect.has_candidate_transcripts(*expected_sequences[sg_id])
            sg_expect.has_transcripts(*expected_sequences[sg_id])

    def test_when_pruning_traverses_two_subgraphs_into_two_transcripts(self, tmpdir):
        # given
        kmer_size = 3
        min_tip_length = 2
        sequences = ['AAATG', 'AAAC', 'ATCCC', 'ATCG']
        seq_recs = [SeqRecord(Seq(seq), id=str(idx),
                              letter_annotations={"phred_quality": [40 for _ in range(len(seq))]})
                    for
                    idx, seq in enumerate(sequences)]
        input_fastq = tmpdir / 'single.fq'
        with open(input_fastq, 'w') as fh:
            SeqIO.write(seq_recs, fh, 'fastq')

        expected_sequences = [['AAATG'], ['ATCCC']]
        expected_kmers = [{'AAA', 'AAT', 'AAC', 'ATG'}, {'ATC', 'GGA', 'CCC', 'CGA'}]
        expected_post_pruning_kmers = [{'AAA', 'AAT', 'ATG'}, {'ATC', 'GGA', 'CCC'}]
        all_expected_kmers = list(itertools.chain(*expected_kmers))
        all_expected_post_pruned_kmers = list(itertools.chain(*expected_post_pruning_kmers))

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
        expect = AbeonaExpectation(out_dir, min_tip_length)
        expect.has_out_graph_with_kmers(*all_expected_kmers)
        expect.has_out_clean_graph_with_kmers(*all_expected_post_pruned_kmers)

        for sg_id in range(2):
            sg_expect = expect.has_subgraph_with_kmers(*expected_post_pruning_kmers[sg_id])
            sg_expect.has_traversal().has_nodes(*expected_post_pruning_kmers[sg_id])
            sg_expect.has_candidate_transcripts(*expected_sequences[sg_id])
            sg_expect.has_transcripts(*expected_sequences[sg_id])

    def test_when_filtering_traverses_three_subgraphs_into_two_transcripts(self, tmpdir):
        # given
        kmer_size = 3
        sequences = [('AAAT', 4), ('ATCC', 4), ('CCCG', 3)]
        sequences = list(itertools.chain(*[[seq for _ in range(num)] for seq, num in sequences]))
        seq_recs = [SeqRecord(Seq(seq), id=str(idx),
                              letter_annotations={"phred_quality": [40 for _ in range(len(seq))]})
                    for
                    idx, seq in enumerate(sequences)]
        input_fastq = tmpdir / 'single.fq'
        with open(input_fastq, 'w') as fh:
            SeqIO.write(seq_recs, fh, 'fastq')

        expected_sequences = [['AAAT'], ['ATCC']]

        all_expected_kmers = {'AAA', 'AAT', 'ATC', 'GGA', 'CCC', 'CCG'}
        expected_clean_kmers = [{'AAA', 'AAT'}, {'ATC', 'GGA'}]
        all_expected_clean_kmers = list(itertools.chain(*expected_clean_kmers))

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
            sg_expect = expect.has_subgraph_with_kmers(*expected_clean_kmers[sg_id])
            sg_expect.has_traversal() \
                .has_nodes(*expected_clean_kmers[sg_id])
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
            sg_expect = expect.has_subgraph_with_kmers(*expected_subgraph_kmers[sg_id])
            sg_expect.has_traversal() \
                .has_nodes(*expected_subgraph_aligned_kmers[sg_id])
            sg_expect.has_candidate_transcripts(*expected_subgraph_sequences[sg_id])
            sg_expect.has_transcripts(*expected_subgraph_sequences[sg_id])
            expect.has_n_subgraphs(1)
