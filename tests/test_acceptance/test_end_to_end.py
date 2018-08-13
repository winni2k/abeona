import gzip
import itertools
from pathlib import Path
from subprocess import check_call

import attr
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from cortexpy.graph.parser.streaming import kmer_string_generator_from_stream, load_cortex_graph
from cortexpy.test.expectation import KmerGraphExpectation
from cortexpy.utils import lexlo


@attr.s(slots=True)
class AbeonaRunner(object):

    def run(self, *args):
        argv = ['python', '-mabeona'] + [str(a) for a in args]
        check_call(argv)

    def assemble(self, *args):
        argv = ['assemble'] + list(args)
        return self.run(*argv)


@attr.s(slots=True)
class FastqBuilder(object):
    out_file = attr.ib()
    records = attr.ib(attr.Factory(list))
    qual = attr.ib(40)

    def with_seq(self, seq_string):
        idx = len(self.records)
        self.records.append(SeqRecord(Seq(seq_string), id=str(idx),
                                      letter_annotations={
                                          "phred_quality": [self.qual for _ in
                                                            range(len(seq_string))]}))
        return self

    def with_seqs(self, *seq_strings):
        for seq_string in seq_strings:
            self.with_seq(seq_string)
        return self

    def build(self):
        with open(self.out_file, 'w') as fh:
            SeqIO.write(self.records, fh, 'fastq')
        return self.out_file


@attr.s(slots=True)
class AbeonaExpectation(object):
    out_dir = attr.ib()
    prune_length = attr.ib(0)
    out_graph = attr.ib(init=False)
    out_clean = attr.ib(init=False)
    subgraphs = attr.ib(init=False)
    traversal_dir = attr.ib(init=False)
    out_clean_tip_pruned = attr.ib(init=False)
    all_transcripts = attr.ib(init=False)

    def __attrs_post_init__(self):
        self.out_dir = Path(self.out_dir)
        self.out_graph = self.out_dir / 'cortex_graph' / 'full.ctx'
        self.out_clean = self.out_graph.with_suffix('.clean.ctx')
        self.out_clean_tip_pruned = self.out_clean.with_suffix(
            f'.min_tip_length_{self.prune_length}.ctx')
        self.traversal_dir = self.out_dir / 'traversals'
        self.all_transcripts = self.out_dir / 'all_transcripts' / 'transcripts.fa.gz'
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
        self._has_kmers_in_graph(kmers, self.out_graph)
        return self

    def has_out_clean_graph_with_kmers(self, *kmers):
        self._has_kmers_in_graph(kmers, self.out_clean)
        return self

    def has_out_tip_pruned_raph_with_kmers(self, *kmers):
        self._has_kmers_in_graph(kmers, self.out_clean_tip_pruned)
        return self

    def has_out_all_transcripts(self, *seqs):
        expected_seqs = [lexlo(s) for s in seqs]
        with gzip.open(self.all_transcripts, 'rt') as fh:
            seqs = list(SeqIO.parse(fh, 'fasta'))
        assert sorted(expected_seqs) == sorted([str(lexlo(rec.seq)) for rec in seqs])
        return self

    def _has_kmers_in_graph(self, kmers, graph):
        with open(graph, 'rb') as fh:
            output_kmers = list(kmer_string_generator_from_stream(fh))
        assert sorted(kmers) == sorted(output_kmers)


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

    def _has_seqs(self, *expected_seq_strings, dir_name, suffix, no_file=False):
        transcripts = self.out_dir / dir_name / f'g{self.sg_id}{suffix}'
        if no_file:
            assert not transcripts.is_file()
            return
        assert transcripts.is_file()
        with gzip.open(str(transcripts), 'rt') as fh:
            seqs = list(SeqIO.parse(fh, 'fasta'))

        for seq in seqs:
            assert seq.id.startswith(f'g{self.sg_id}_')

        seq_strings = {str(s.seq) for s in seqs} | {str(s.seq.reverse_complement()) for s in seqs}
        for expected_seq_string in expected_seq_strings:
            assert expected_seq_string in seq_strings
        assert len(set(expected_seq_strings)) == len(seqs)

    def has_candidate_transcripts(self, *seq_strings, no_file=False):
        self._has_seqs(*seq_strings, dir_name='candidate_transcripts',
                       suffix='.candidate_transcripts.fa.gz', no_file=no_file)
        return self

    def has_no_candidate_transcripts(self):
        self.has_candidate_transcripts()
        return self

    def has_transcripts(self, *seq_strings, no_file=False):
        self._has_seqs(*seq_strings, dir_name='transcripts', suffix='.transcripts.fa.gz',
                       no_file=no_file)

    def has_no_transcripts(self):
        self.has_transcripts(no_file=True)
        return self


class TestAssemble(object):
    def test_raises_without_input(self, tmpdir):
        b = FastqBuilder(tmpdir / 'input.fastq')
        input_fastq = b.build()
        out_dir = Path(tmpdir) / 'abeona'
        args = [
            '--fastx-single', input_fastq,
            '--kallisto-fastx-single', input_fastq,
            '--kallisto-fragment-length', 31,
            '--kallisto-sd', 0.1,
            '--bootstrap-samples', 100,
            '--out-dir', out_dir,
            '--kmer-size', 47,
            '--min-unitig-coverage', 0,
        ]

        # when/then
        with pytest.raises(Exception):
            AbeonaRunner().assemble(*args)

    def test_traverses_two_subgraphs_into_two_transcripts(self, tmpdir):
        # given
        kmer_size = 3
        b = FastqBuilder(tmpdir / 'single.fq')
        b.with_seqs('AAAT', 'ATCC')
        input_fastq = b.build()

        expected_sequences = [['AAAT'], ['ATCC']]
        expected_kmers = [{'AAA', 'AAT'}, {'ATC', 'GGA'}]
        all_expected_kmers = list(itertools.chain(*expected_kmers))

        out_dir = Path(tmpdir) / 'abeona'
        args = ['--fastx-single', input_fastq,
                '--kallisto-fastx-single', input_fastq,
                '--kallisto-fragment-length', 3,
                '--kallisto-sd', 0.1,
                '--bootstrap-samples', 100,
                '--out-dir', out_dir,
                '--kmer-size', kmer_size,
                '--min-unitig-coverage', 0, ]

        # when
        AbeonaRunner().assemble(*args)

        # then
        expect = AbeonaExpectation(out_dir)
        expect.has_out_graph_with_kmers(*all_expected_kmers)
        expect.has_out_clean_graph_with_kmers(*all_expected_kmers)

        for sg_id in range(2):
            sg_expect = expect.has_subgraph_with_kmers(*expected_kmers[sg_id])
            sg_expect.has_traversal().has_nodes(*expected_kmers[sg_id])
            sg_expect.has_transcripts(*expected_sequences[sg_id])

        expect.has_out_all_transcripts('AAAT', 'ATCC')

    def test_traverses_two_subgraphs_of_two_transcripts_into_four_transcripts(self, tmpdir):
        # given
        kmer_size = 3
        b = FastqBuilder(tmpdir / 'single.fq')
        for _ in range(4):
            b.with_seqs('AAAT', 'AAAC', 'ATCC', 'ATCA')

        input_fastq = b.build()

        out_dir = Path(tmpdir) / 'abeona'

        # when
        AbeonaRunner().assemble('--fastx-single', input_fastq,
                                '--kallisto-fastx-single', input_fastq,
                                '--kallisto-fragment-length', 3,
                                '--kallisto-sd', 0.1,
                                '--bootstrap-samples', 100,
                                '--out-dir', out_dir,
                                '--kmer-size', kmer_size,
                                '--min-unitig-coverage', 0)

        # then
        expect = AbeonaExpectation(out_dir)
        expect.has_out_graph_with_kmers('AAA', 'AAT', 'AAC', 'ATC', 'GGA', 'TCA')
        expect.has_out_clean_graph_with_kmers('AAA', 'AAT', 'AAC', 'ATC', 'GGA', 'TCA')

        sg = expect.has_subgraph_with_kmers('AAA', 'AAT', 'AAC')
        sg.has_traversal().has_nodes('AAA', 'AAT', 'AAC')
        sg.has_candidate_transcripts('AAAT', 'AAAC')
        sg.has_transcripts('AAAT', 'AAAC')

        sg = expect.has_subgraph_with_kmers('ATC', 'GGA', 'TCA')
        sg.has_traversal().has_nodes('ATC', 'GGA', 'TCA')
        sg.has_candidate_transcripts('ATCC', 'ATCA')
        sg.has_transcripts('ATCC', 'ATCA')

    def test_when_pruning_traverses_two_subgraphs_into_two_transcripts(self, tmpdir):
        # given
        min_tip_length = 2
        b = FastqBuilder(tmpdir / 'single.fq')
        for seq in ['AAATG', 'AAAC', 'ATCCC', 'ATCG']:
            b.with_seq(seq)
        input_fastq = b.build()

        expected_sequences = [['AAATG'], ['ATCCC']]
        expected_kmers = [{'AAA', 'AAT', 'AAC', 'ATG'}, {'ATC', 'GGA', 'CCC', 'CGA'}]
        expected_post_pruning_kmers = [{'AAA', 'AAT', 'ATG'}, {'ATC', 'GGA', 'CCC'}]
        all_expected_kmers = list(itertools.chain(*expected_kmers))
        all_expected_post_pruned_kmers = list(itertools.chain(*expected_post_pruning_kmers))

        out_dir = Path(tmpdir) / 'abeona'
        args = ['--fastx-single', input_fastq,
                '--kallisto-fastx-single', input_fastq,
                '--kallisto-fragment-length', 3,
                '--kallisto-sd', 0.1,
                '--bootstrap-samples', 100,
                '--out-dir', out_dir,
                '--kmer-size', 3,
                '--min-tip-length', min_tip_length,
                '--min-unitig-coverage', 0,
                '--quiet']

        # when
        AbeonaRunner().assemble(*args)

        # then
        expect = AbeonaExpectation(out_dir, min_tip_length)
        expect.has_out_graph_with_kmers(*all_expected_kmers)
        expect.has_out_clean_graph_with_kmers(*all_expected_kmers)
        expect.has_out_tip_pruned_raph_with_kmers(*all_expected_post_pruned_kmers)

        for sg_id in range(2):
            sg_expect = expect.has_subgraph_with_kmers(*expected_post_pruning_kmers[sg_id])
            sg_expect.has_traversal().has_nodes(*expected_post_pruning_kmers[sg_id])
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
        args = ['--fastx-single', input_fastq,
                '--kallisto-fastx-single', input_fastq,
                '--kallisto-fragment-length', 3,
                '--kallisto-sd', 0.1,
                '--bootstrap-samples', 100,
                '--out-dir', out_dir,
                '--kmer-size', kmer_size]

        # when
        AbeonaRunner().assemble(*args)

        # then
        expect = AbeonaExpectation(out_dir)
        expect.has_out_graph_with_kmers(*all_expected_kmers)
        expect.has_out_clean_graph_with_kmers(*all_expected_clean_kmers)

        for sg_id in range(2):
            sg_expect = expect.has_subgraph_with_kmers(*expected_clean_kmers[sg_id])
            sg_expect.has_traversal() \
                .has_nodes(*expected_clean_kmers[sg_id])
            sg_expect.has_transcripts(*expected_sequences[sg_id])

    @pytest.mark.skip(
        reason='Graphs with single unitigs are not passed through Kallisto at this time')
    def test_ignores_polyA_tails_with_kmer_size_47(self, tmpdir):
        # given
        b = FastqBuilder(tmpdir / 'input.fastq')

        for _ in range(4):
            b.with_seq('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')
            b.with_seq('ATATATATATATATATATATATATATATATATATATATATATATATA')

        input_fastq = b.build()
        out_dir = Path(tmpdir) / 'abeona'
        args = [
            '--fastx-single', input_fastq,
            '--kallisto-fastx-single', input_fastq,
            '--kallisto-fragment-length', 31,
            '--kallisto-sd', 0.1,
            '--bootstrap-samples', 100,
            '--out-dir', out_dir,
            '--kmer-size', 47,
            '--min-unitig-coverage', 0,
        ]

        # when
        AbeonaRunner().assemble(*args)

        # then
        expect = AbeonaExpectation(out_dir)
        expect.has_out_graph_with_kmers('ATATATATATATATATATATATATATATATATATATATATATATATA',
                                        'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')

        expect \
            .has_subgraph_with_kmers('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA') \
            .has_no_transcripts()
        expect \
            .has_subgraph_with_kmers('ATATATATATATATATATATATATATATATATATATATATATATATA') \
            .has_transcripts('ATATATATATATATATATATATATATATATATATATATATATATATA')
        expect.has_n_subgraphs(2)


class TestMaxPaths(object):
    def test_when_using_max_paths_traverses_two_subgraphs_into_one_transcript(self, tmpdir):
        # given
        kmer_size = 3

        sequences = ['AAAT', 'AAAC', 'ATCC']
        b = FastqBuilder(tmpdir / 'input.fastq')
        [b.with_seq(seq) for seq in sequences]
        input_fastq = b.build()

        out_dir = Path(tmpdir) / 'abeona'
        args = ['--fastx-single', input_fastq,
                '--kallisto-fastx-single', input_fastq,
                '--kallisto-fragment-length', 3,
                '--kallisto-sd', 0.1,
                '--bootstrap-samples', 100,
                '--out-dir', out_dir,
                '--kmer-size', kmer_size,
                '--min-unitig-coverage', 0,
                '--max-paths-per-subgraph', 1]

        # when
        AbeonaRunner().assemble(*args)

        # then
        expect = AbeonaExpectation(out_dir)
        expect.has_out_graph_with_kmers('ATC', 'GGA', 'AAA', 'AAT', 'AAC')
        expect.has_out_clean_graph_with_kmers('ATC', 'GGA', 'AAA', 'AAT', 'AAC')

        # correct graph
        sg_expect = expect.has_subgraph_with_kmers('ATC', 'GGA')
        sg_expect.has_traversal().has_nodes('ATC', 'GGA')
        sg_expect.has_transcripts('ATCC')

        # graph with too many paths
        expect = expect.has_subgraph_with_kmers('AAA', 'AAC', 'AAT')
        expect.has_no_transcripts()


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

        b = FastqBuilder(tmpdir / 'single.fq')
        for seq in sequences:
            b.with_seq(seq)
        input_fastq = b.build()

        expected_kmers = [{'AAA', 'AAT'}, {'ATC', 'GGA'}]
        all_expected_kmers = list(itertools.chain(*expected_kmers))

        expected_subgraph_sequences = [['AAAT']]
        expected_subgraph_kmers = [{'AAA', 'AAT'}]
        expected_subgraph_aligned_kmers = [{'AAA', 'AAT'}, {'ATC', 'TCC'}]

        out_dir = Path(tmpdir) / 'abeona'
        args = ['--initial-contigs', inital_contigs_fasta,
                '--fastx-single', input_fastq,
                '--kallisto-fastx-single', input_fastq,
                '--kallisto-fragment-length', 3,
                '--kallisto-sd', 0.1,
                '--bootstrap-samples', 100,
                '--out-dir', out_dir,
                '--kmer-size', kmer_size,
                '--min-unitig-coverage', 0, ]

        # when
        AbeonaRunner().assemble(*args)

        # then
        expect = AbeonaExpectation(out_dir)
        expect.has_out_graph_with_kmers(*all_expected_kmers)
        expect.has_out_clean_graph_with_kmers(*all_expected_kmers)

        for sg_id in range(1):
            sg_expect = expect.has_subgraph_with_kmers(*expected_subgraph_kmers[sg_id])
            sg_expect.has_traversal() \
                .has_nodes(*expected_subgraph_aligned_kmers[sg_id])
            sg_expect.has_transcripts(*expected_subgraph_sequences[sg_id])
            expect.has_n_subgraphs(1)
