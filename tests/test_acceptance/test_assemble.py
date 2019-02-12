import gzip
import itertools
import re
import zipfile
from pathlib import Path
from subprocess import check_call

import attr
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from cortexpy.graph.parser.streaming import kmer_string_generator_from_stream, load_cortex_graph
from cortexpy.links import Links
from cortexpy.test.expectation import KmerGraphExpectation
from cortexpy.utils import lexlo


@attr.s(slots=True)
class AbeonaRunner:

    def run(self, *args):
        argv = ['abeona'] + [str(a) for a in args]
        check_call(argv)

    def assemble(self, *args, no_cleanup=True):
        argv = ['assemble'] + list(args)
        if '--min-tip-length' not in argv:
            argv += ['--min-tip-length', '0']
        if no_cleanup:
            argv.append('--no-cleanup')
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
class PairedFastqBuilder(object):
    out_file1 = attr.ib()
    out_file2 = attr.ib()
    records1 = attr.ib(attr.Factory(list))
    records2 = attr.ib(attr.Factory(list))
    qual = attr.ib(40)

    def with_pair(self, forward, reverse):
        idx = len(self.records1)
        quals_f = [self.qual for _ in range(len(forward))]
        quals_r = [self.qual for _ in range(len(reverse))]
        self.records1.append(
            SeqRecord(Seq(forward), id=f'{idx}', letter_annotations={"phred_quality": quals_f})
        )
        self.records2.append(
            SeqRecord(Seq(reverse), id=f'{idx}', letter_annotations={"phred_quality": quals_r})
        )
        return self

    def build(self):
        with open(self.out_file1, 'w') as fh:
            SeqIO.write(self.records1, fh, 'fastq')
        with open(self.out_file2, 'w') as fh:
            SeqIO.write(self.records2, fh, 'fastq')
        return self.out_file1, self.out_file2


@attr.s(slots=True)
class AbeonaExpectation(object):
    out_dir = attr.ib()
    prune_length = attr.ib(0)
    out_graph = attr.ib(init=False)
    out_clean = attr.ib(init=False)
    subgraphs = attr.ib(init=False)
    subgraphs_zip = attr.ib(init=False)
    traversal_dir = attr.ib(init=False)
    out_clean_tip_pruned = attr.ib(init=False)
    all_transcripts = attr.ib(init=False)
    skipped_subgraphs = attr.ib(init=False)

    def __attrs_post_init__(self):
        self.out_dir = Path(self.out_dir)
        self.out_graph = self.out_dir / 'cortex_graph' / 'full.ctx'
        self.out_clean = self.out_graph.with_suffix('.clean.ctx')
        self.out_clean_tip_pruned = self.out_clean.with_suffix(
            f'.min_tip_length_{self.prune_length}.ctx')
        self.traversal_dir = self.out_dir / 'traversals'
        self.all_transcripts = self.out_dir / 'transcripts.fa'
        self.skipped_subgraphs = self.out_dir / 'skipped_subgraphs' / 'skipped_subgraphs.txt'

        self.subgraphs_zip = self.traversal_dir / 'subgraphs.zip'
        self.subgraphs = []
        sg_pattern = re.compile('g(\d+)')
        try:
            with zipfile.ZipFile(self.subgraphs_zip, 'r') as zip_fh:
                for sg in zip_fh.namelist():
                    sg_id = sg_pattern.search(sg).group(1)
                    self.subgraphs.append((sg_id, sg))
        except FileNotFoundError:
            pass

    def has_subgraph(self, sg_id):
        return AbeonaSubgraphExpectation.from_abeona_expectation(self, sg_id)

    def has_subgraph_with_kmers(self, *kmers):
        kmers = set(kmers)
        assert len(kmers) > 0
        for sg_id, sg in self.subgraphs:
            with zipfile.ZipFile(self.subgraphs_zip, 'r') as zfh:
                with zfh.open(sg, 'r') as fh:
                    subgraph_kmers = set(lexlo(n) for n in load_cortex_graph(fh))
            if len(kmers & subgraph_kmers) > 0:
                assert kmers == subgraph_kmers
                return AbeonaSubgraphExpectation.from_abeona_expectation(self, sg_id)
        assert kmers == set()
        return self

    def has_subgraph_with_kmer(self, kmer):
        for sg_id, sg in self.subgraphs:
            with zipfile.ZipFile(self.subgraphs_zip, 'r') as zfh:
                with zfh.open(sg, 'r') as fh:
                    subgraph_kmers = set(lexlo(n) for n in load_cortex_graph(fh))
            if kmer in subgraph_kmers:
                return AbeonaSubgraphExpectation.from_abeona_expectation(self, sg_id)
        assert False, f'Could not find subgraph with kmer: {kmer}'

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
        if not self.all_transcripts.is_file():
            assert 0 == len(seqs)
            return self
        expected_seqs = [lexlo(s) for s in seqs]
        with open(self.all_transcripts, 'rt') as fh:
            seqs = [str(lexlo(rec.seq)) for rec in SeqIO.parse(fh, 'fasta')]
        assert sorted(expected_seqs) == sorted(seqs)
        return self

    def has_out_all_transcript_description_matching(self, regex):
        with open(self.all_transcripts, 'rt') as fh:
            descs = [rec.description for rec in SeqIO.parse(fh, 'fasta')]
        desc_re = re.compile(regex)
        assert any(desc_re.search(d) for d in descs)
        return self

    def has_n_missing_subgraphs(self, n):
        assert n == len(self.skipped_subgraphs.read_text().splitlines())
        return self

    def _has_kmers_in_graph(self, kmers, graph):
        with open(graph, 'rb') as fh:
            output_kmers = list(kmer_string_generator_from_stream(fh))
        assert sorted(kmers) == sorted(output_kmers)

    def has_out_dirs(self, *dirs):
        assert sorted(dirs) == sorted(d.name for d in self.out_dir.iterdir() if d.is_dir())
        return self


@attr.s(slots=True)
class AbeonaSubgraphExpectation:
    out_dir = attr.ib()
    sg_id = attr.ib()
    sg_traversal_name = attr.ib()
    subgraphs_zip = attr.ib()

    def has_traversal_graph_with_kmers(self, *kmers):
        with zipfile.ZipFile(self.subgraphs_zip, 'r') as zfh:
            with zfh.open(self.sg_traversal_name, 'r') as fh:
                output_kmers = set(kmer_string_generator_from_stream(fh))
        assert set(kmers) == output_kmers
        return self

    def has_traversal(self):
        with zipfile.ZipFile(self.subgraphs_zip, 'r') as zfh:
            with zfh.open(self.sg_traversal_name, 'r') as fh:
                return KmerGraphExpectation(load_cortex_graph(fh))

    def _has_transcripts(self, *expected_seq_strings, dir_name, suffix, id_check=True):
        transcripts = self.out_dir / dir_name / f'g{self.sg_id}{suffix}'
        if len(expected_seq_strings) == 0:
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
        assert len(expected_seq_strings) == len(seqs)

    def _has_reads(self, *expected_seq_strings, dir_name, suffix, file_type='fastq'):
        reads = self.out_dir / dir_name / f'g{self.sg_id}{suffix}'
        if len(expected_seq_strings) == 0:
            assert not reads.is_file()
            return
        assert reads.is_file()
        with open(str(reads), 'rt') as fh:
            seqs = list(SeqIO.parse(fh, file_type))

        seq_strings = [str(s.seq) for s in seqs]
        assert sorted(expected_seq_strings) == sorted(seq_strings)

    def has_candidate_transcripts(self, *seq_strings):
        self._has_transcripts(*seq_strings, dir_name='candidate_transcripts',
                              suffix='.candidate_transcripts.fa.gz')
        return self

    def has_reads_assigned(self, reads, reads2=None):
        self._has_reads(*reads, dir_name='reads_assigned_to_subgraphs', suffix='.1.fa',
                        file_type='fasta')
        if reads2 is not None:
            self._has_reads(*reads2, dir_name='reads_assigned_to_subgraphs', suffix='.2.fa',
                            file_type='fasta')
        return self

    def has_no_candidate_transcripts(self):
        self.has_candidate_transcripts()
        return self

    def has_transcripts(self, *seq_strings):
        self._has_transcripts(*seq_strings, dir_name='transcripts', suffix='.transcripts.fa.gz')

    def has_no_transcripts(self):
        self.has_transcripts()
        return self

    def has_links_for_kmers(self, *kmers):
        with gzip.open(self.out_dir / 'links' / f'g{self.sg_id}.ctp.gz', 'rb') as fh:
            links = Links.from_binary_stream(fh)
            for kmer in kmers:
                assert lexlo(kmer) == kmer
                assert kmer in links.body.keys()

    @classmethod
    def from_abeona_expectation(cls, abeona_expectation, sg_id):
        return cls(abeona_expectation.out_dir, sg_id,
                   subgraphs_zip=abeona_expectation.subgraphs_zip,
                   sg_traversal_name=f'g{sg_id}.traverse.ctx')


class TestAssemble:
    def test_raises_without_input(self, tmpdir):
        b = FastqBuilder(tmpdir / 'input.fastq')
        input_fastq = b.build()
        out_dir = Path(tmpdir) / 'abeona'
        args = [
            '--fastx-single', input_fastq,
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

    def test_cleanup(self, tmpdir):
        b = FastqBuilder(tmpdir / 'single.fq')
        b.with_seq('AAAC')

        out_dir = Path(tmpdir) / 'abeona'
        args = ['--fastx-single', b.build(),
                '--kallisto-fragment-length', 3,
                '--kallisto-sd', 0.1,
                '--bootstrap-samples', 10,
                '--out-dir', out_dir,
                '--kmer-size', 3,
                '--min-unitig-coverage', 0, ]

        # when
        AbeonaRunner().assemble(*args, no_cleanup=False)

        # then
        expect = AbeonaExpectation(out_dir)
        expect.has_out_all_transcripts('AAAC')
        expect.has_out_dirs('.nextflow')

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
                '--kallisto-fragment-length', 3,
                '--kallisto-sd', 0.1,
                '--bootstrap-samples', 10,
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
        b = FastqBuilder(tmpdir / 'single.fq')
        for _ in range(4):
            b.with_seqs('AAAAAT',
                        'AAAATT',
                        'AAAAAC',
                        'AAAACA',
                        'AAATCC',
                        'AATCCC',
                        'AAATCT',
                        'AATCTT')

        input_fastq = b.build()

        out_dir = Path(tmpdir) / 'abeona'
        args = ['--fastx-single', input_fastq,
                '--kallisto-fragment-length', 3,
                '--kallisto-sd', 0.1,
                '--bootstrap-samples', 100,
                '--out-dir', out_dir,
                '--kmer-size', 5,
                '--min-unitig-coverage', 0,
                '--record-buffer-size', 1]

        # when
        AbeonaRunner().assemble(*args)

        # then
        expect = AbeonaExpectation(out_dir)
        all_kmers = [
            'AAAAA',
            'AAAAT',
            'AAATT',
            'AAAAC',
            'AAACA',
            'AAATC',
            'AATCC',
            'ATCCC',
            'AATCT',
            'AAGAT',
        ]
        expect.has_out_graph_with_kmers(*all_kmers)
        expect.has_out_clean_graph_with_kmers(*all_kmers)

        sg = expect.has_subgraph_with_kmers('AAAAA', 'AAAAT', 'AAATT', 'AAAAC', 'AAACA')
        sg.has_traversal().has_nodes('AAAAA', 'AAAAT', 'AAATT', 'AAAAC', 'AAACA')
        sg.has_candidate_transcripts('AAAAATT', 'AAAAACA')
        sg.has_reads_assigned(['AAAAAT', 'AAAATT', 'AAAAAC', 'AAAACA'] * 4)
        sg.has_transcripts('AAAAATT', 'AAAAACA')

        sg = expect.has_subgraph_with_kmers('AAATC', 'AATCC', 'ATCCC', 'AATCT', 'AAGAT')
        sg.has_traversal().has_nodes('AAATC', 'AATCC', 'ATCCC', 'AATCT', 'AAGAT')
        sg.has_candidate_transcripts('AAATCCC', 'AAATCTT')
        sg.has_transcripts('AAATCCC', 'AAATCTT')
        sg.has_reads_assigned(['AAATCC', 'AATCCC', 'AAATCT', 'AATCTT'] * 4)

        expect.has_out_all_transcripts('AAAAATT', 'AAAAACA', 'AAATCCC', 'AAATCTT')

    def test_ignores_subgraph_with_unitig_shorter_than_reads(self, tmpdir):

        # given
        b = FastqBuilder(tmpdir / 'single.fq')
        for _ in range(1):
            b.with_seqs(
                'TAAAAAT',
                'TAAAAATT',
                'CAAAAAT',
                'CAAAAATC',
            )

        input_fastq = b.build()

        out_dir = Path(tmpdir) / 'abeona'
        args = ['--fastx-single', input_fastq,
                '--kallisto-fragment-length', 7,
                '--kallisto-sd', 0.1,
                '--bootstrap-samples', 100,
                '--out-dir', out_dir,
                '--kmer-size', 5,
                '--min-unitig-coverage', 2]

        # when
        AbeonaRunner().assemble(*args)

        # then
        expect = AbeonaExpectation(out_dir)
        expect \
            .has_subgraph(0) \
            .has_traversal_graph_with_kmers('AAAAA', 'AAAAT', 'TAAAA', 'CAAAA') \
            .has_reads_assigned(['TAAAAATT', 'CAAAAATC', 'TAAAAAT', 'CAAAAAT']) \
            .has_no_candidate_transcripts()

    def test_traverses_tangle_into_two_candidate_transcripts(self, tmpdir):

        # given
        b = FastqBuilder(tmpdir / 'single.fq')
        for _ in range(4):
            b.with_seqs('TAAAAAT',
                        'AAAATA',
                        'GAAAAAC',
                        'AAAACA')

        input_fastq = b.build()

        out_dir = Path(tmpdir) / 'abeona'
        args = ['--fastx-single', input_fastq,
                '--kallisto-fragment-length', 7,
                '--kallisto-sd', 0.1,
                '--bootstrap-samples', 100,
                '--out-dir', out_dir,
                '--kmer-size', 5,
                '--min-unitig-coverage', 0]

        # when
        AbeonaRunner().assemble(*args)

        # then
        expect = AbeonaExpectation(out_dir)
        expect \
            .has_subgraph(0) \
            .has_candidate_transcripts('TAAAAATA',
                                       'GAAAAACA')
        expect.has_out_all_transcripts('TAAAAATA',
                                       'GAAAAACA')

    def test_ignores_subgraph_with_two_junctions(self, tmpdir):

        # given
        b = FastqBuilder(tmpdir / 'single.fq')
        for _ in range(4):
            b.with_seqs('TAAAAAT',
                        'AAAATA',
                        'GAAAAAC',
                        'AAAACA')

        input_fastq = b.build()

        out_dir = Path(tmpdir) / 'abeona'
        args = ['--fastx-single', input_fastq,
                '--kallisto-fragment-length', 7,
                '--kallisto-sd', 0.1,
                '--bootstrap-samples', 100,
                '--out-dir', out_dir,
                '--kmer-size', 5,
                '--max-junctions', 1,
                '--min-unitig-coverage', 0]

        # when
        AbeonaRunner().assemble(*args)

        # then
        expect = AbeonaExpectation(out_dir)
        expect.has_n_missing_subgraphs(1)

    def test_with_read_pairs_traverses_graph_of_two_transcripts_into_two_transcripts(self, tmpdir):

        # given
        kmer_size = 29
        b = PairedFastqBuilder(tmpdir / 'forward.fq', tmpdir / 'reverse.fq')
        read_pairs = [
            ['AAAAAAAAAAAAAAAAAAAAAAAAAACTCAAAAAAAAAAAAAAAAAAAAAACCCC',
             'CTCAAAAAAAAAAAAAAAAAAAAAACCCCAAAAAAAAAAAAAAAAAAAAAAAAAAA'],
            ['AAAAAAAAAAAAAAAAAAAAAAAAAACACAAAAAAAAAAAAAAAAAAAAAACCCC',
             'CACAAAAAAAAAAAAAAAAAAAAAACCCCAAAAAAAAAAAAAAAAAAAAAAAAAAA']
        ]
        for _ in range(4):
            b.with_pair(*read_pairs[0])
            b.with_pair(*read_pairs[1])

        expected_transcripts = [
            'AAAAAAAAAAAAAAAAAAAAAAAAAACTCAAAAAAAAAAAAAAAAAAAAAACCCCAAAAAAAAAAAAAAAAAAAAAAAAAAA',
            'AAAAAAAAAAAAAAAAAAAAAAAAAACACAAAAAAAAAAAAAAAAAAAAAACCCCAAAAAAAAAAAAAAAAAAAAAAAAAAA'
        ]

        forward, reverse = b.build()

        out_dir = Path(tmpdir) / 'abeona'
        args = ['--fastx-forward', forward,
                '--fastx-reverse', reverse,
                '--bootstrap-samples', 100,
                '--out-dir', out_dir,
                '--kmer-size', kmer_size,
                '--min-unitig-coverage', 0]

        # when
        AbeonaRunner().assemble(*args)

        # then
        expect = AbeonaExpectation(out_dir)
        sg = expect.has_subgraph(0)
        sg.has_candidate_transcripts(*expected_transcripts)
        sg.has_reads_assigned(['AAAAAAAAAAAAAAAAAAAAAAAAAACTCAAAAAAAAAAAAAAAAAAAAAACCCC',
                               'AAAAAAAAAAAAAAAAAAAAAAAAAACACAAAAAAAAAAAAAAAAAAAAAACCCC'] * 4,
                              ['CTCAAAAAAAAAAAAAAAAAAAAAACCCCAAAAAAAAAAAAAAAAAAAAAAAAAAA',
                               'CACAAAAAAAAAAAAAAAAAAAAAACCCCAAAAAAAAAAAAAAAAAAAAAAAAAAA'] * 4)
        sg.has_transcripts(*expected_transcripts)

        expect.has_out_all_transcripts(*expected_transcripts)

    def test_with_read_pairs_traverses_graph_of_two_transcripts_into_two_graphs(self, tmpdir):

        # given
        kmer_size = 5
        b = PairedFastqBuilder(tmpdir / 'forward.fq', tmpdir / 'reverse.fq')
        read_pairs = [
            ['AAAAAC', 'AAAACC'],
            ['AAAAAC', 'AAAACT'],
            ['CCCCCT', 'CCCCTT'],
            ['CCCCCT', 'CCCCTA'],
        ]
        for _ in range(4):
            for rp in read_pairs:
                b.with_pair(*rp)

        forward, reverse = b.build()

        out_dir = Path(tmpdir) / 'abeona'
        args = ['--fastx-forward', forward,
                '--fastx-reverse', reverse,
                '--bootstrap-samples', 100,
                '--out-dir', out_dir,
                '--kmer-size', kmer_size,
                '--min-unitig-coverage', 0]

        # when
        AbeonaRunner().assemble(*args)

        # then
        expect = AbeonaExpectation(out_dir)
        sg = expect.has_subgraph_with_kmers('AAAAA', 'AAAAC', 'AAACC', 'AAACT')
        sg.has_reads_assigned(['AAAAAC'] * 8, ['AAAACC', 'AAAACT'] * 4)

        sg = expect.has_subgraph_with_kmers('CCCCC', 'AGGGG', 'AAGGG', 'CCCTA')
        sg.has_reads_assigned(['CCCCCT'] * 8, ['CCCCTT', 'CCCCTA'] * 4)

        expect.has_out_all_transcripts('AAAAACC', 'AAAAACT', 'AAGGGGG', 'CCCCCTA')

    @pytest.mark.parametrize('prune_tips_with_mccortex', [True, False])
    def test_when_pruning_traverses_two_subgraphs_into_two_transcripts(self, tmpdir,
                                                                       prune_tips_with_mccortex):
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
                '--kallisto-fragment-length', 3,
                '--kallisto-sd', 0.1,
                '--bootstrap-samples', 100,
                '--out-dir', out_dir,
                '--kmer-size', 3,
                '--min-tip-length', min_tip_length,
                '--min-unitig-coverage', 0,
                '--quiet']
        if prune_tips_with_mccortex:
            args += ['--prune-tips-with-mccortex']

        # when
        AbeonaRunner().assemble(*args)

        # then
        expect = AbeonaExpectation(out_dir, min_tip_length)
        expect.has_out_graph_with_kmers(*all_expected_kmers)
        if prune_tips_with_mccortex:
            expect.has_out_clean_graph_with_kmers(*all_expected_post_pruned_kmers)
        else:
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


class TestLinks:
    @pytest.mark.parametrize('no_links', [True, False])
    def test_creates_two_sets_of_links_for_two_tangles_in_two_subgraphs(self, tmpdir, no_links):
        # given
        b = FastqBuilder(tmpdir / 'single.fq')
        for _ in range(4):
            b.with_seqs('TAAAAAT',
                        'AAAATA',
                        'GAAAAAC',
                        'AAAACA')
            b.with_seqs('TCCCCCT',
                        'CCCCTT',
                        'GCCCCCA',
                        'CCCCAA')

        input_fastq = b.build()

        out_dir = Path(tmpdir) / 'abeona'
        args = [
            '--fastx-single', input_fastq,
            '--kallisto-fragment-length', 7,
            '--kallisto-sd', 0.1,
            '--bootstrap-samples', 100,
            '--out-dir', out_dir,
            '--kmer-size', 5,
            '--min-unitig-coverage', 0,
        ]
        if no_links:
            args.append('--no-links')

        # when
        AbeonaRunner().assemble(*args)

        # then
        expect = AbeonaExpectation(out_dir)

        sg = expect.has_subgraph_with_kmer('AAAAA')
        sg.has_links_for_kmers('TAAAA', 'GAAAA', 'AAAAT', 'AAAAC')
        if no_links:
            sg.has_candidate_transcripts('TAAAAATA', 'GAAAAACA', 'GAAAAATA', 'TAAAAACA')
        else:
            sg.has_candidate_transcripts('TAAAAATA', 'GAAAAACA')

        sg = expect.has_subgraph_with_kmer('CCCCC')
        sg.has_links_for_kmers('GGGGA', 'GCCCC', 'AGGGG', 'CCCCA')

        if no_links:
            sg.has_candidate_transcripts('TCCCCCTT', 'GCCCCCAA', 'GCCCCCTT', 'TCCCCCAA')
        else:
            sg.has_candidate_transcripts('TCCCCCTT', 'GCCCCCAA')

        expect.has_out_all_transcripts('TAAAAATA', 'GAAAAACA', 'TCCCCCTT', 'GCCCCCAA')


class TestMaxPaths:
    def test_when_using_max_paths_traverses_two_subgraphs_into_one_transcript(self, tmpdir):
        # given
        kmer_size = 3

        sequences = ['AAAT', 'AAAC', 'ATCC']
        b = FastqBuilder(tmpdir / 'input.fastq')
        [b.with_seq(seq) for seq in sequences]
        input_fastq = b.build()

        out_dir = Path(tmpdir) / 'abeona'
        args = ['--fastx-single', input_fastq,
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
        expect.has_out_all_transcripts('ATCC')

        # correct graph
        sg_expect = expect.has_subgraph_with_kmers('ATC', 'GGA')
        sg_expect.has_traversal().has_nodes('ATC', 'GGA')
        sg_expect.has_transcripts('ATCC')

        # graph with too many paths
        skipped_expect = expect.has_subgraph_with_kmers('AAA', 'AAC', 'AAT')
        skipped_expect.has_no_transcripts()
        expect.has_n_missing_subgraphs(1)

    def test_traverses_cycle_into_no_transcripts(self, tmpdir):
        # given
        kmer_size = 3
        b = FastqBuilder(tmpdir / 'single.fq')
        b.with_seqs('AAACAAA')
        input_fastq = b.build()

        out_dir = Path(tmpdir) / 'abeona'
        args = ['--fastx-single', input_fastq,
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
        expect.has_out_graph_with_kmers('AAA', 'AAC', 'ACA', 'CAA')
        expect.has_out_clean_graph_with_kmers('AAA', 'AAC', 'ACA', 'CAA')

        expect.has_out_all_transcripts()


class TestInitialSeqsFasta:
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
            b.with_seq(seq)
        input_fastq = b.build()

        expected_kmers = [{'AAA', 'AAT'}, {'ATC', 'GGA'}]
        all_expected_kmers = list(itertools.chain(*expected_kmers))

        out_dir = Path(tmpdir) / 'abeona'
        args = ['--initial-contigs', inital_contigs_fasta,
                '--fastx-single', input_fastq,
                '--kallisto-fragment-length', 3,
                '--kallisto-sd', 0.1,
                '--bootstrap-samples', 100,
                '--out-dir', out_dir,
                '--kmer-size', kmer_size,
                '--min-unitig-coverage', 0,
                ]

        # when
        AbeonaRunner().assemble(*args)

        # then
        expect = AbeonaExpectation(out_dir)
        expect.has_out_graph_with_kmers(*all_expected_kmers)
        expect.has_out_clean_graph_with_kmers(*all_expected_kmers)

        sg_expect = expect.has_subgraph_with_kmers('AAA', 'AAT')
        sg_expect.has_traversal() \
            .has_nodes('AAA', 'AAT')
        sg_expect.has_transcripts('AAAT')

        expect.has_n_subgraphs(1)
        expect.has_out_all_transcripts('AAAT')

    def test_with_extra_start_kmer_assembles_circular_graph_into_single_sequence(self, tmpdir):
        # given
        kmer_size = 3
        sequences = ['CCCACCC']

        inital_contigs_fasta = tmpdir / 'initial-contigs.fa'
        SeqIO.write([SeqRecord(Seq('CCC'), id='0')], str(inital_contigs_fasta), 'fasta')

        b = FastqBuilder(tmpdir / 'single.fq')
        for seq in sequences:
            b.with_seq(seq)
            b.with_seq(seq)
        input_fastq = b.build()

        out_dir = Path(tmpdir) / 'abeona'
        args = ['--initial-contigs', inital_contigs_fasta,
                '--fastx-single', input_fastq,
                '--kallisto-fragment-length', len(sequences[0]),
                '--kallisto-sd', 0.1,
                '--bootstrap-samples', 100,
                '--out-dir', out_dir,
                '--kmer-size', kmer_size,
                '--min-unitig-coverage', 0,
                '--extra-start-kmer', 'CAC'
                ]

        # when
        AbeonaRunner().assemble(*args)

        # then
        expect = AbeonaExpectation(out_dir)
        expect.has_out_all_transcripts('CACCCA')


class TestParametrizeFiltering:
    @pytest.mark.parametrize('estimated_count_threshold', [0, 5])
    def test_with_bootstrap_count_threshold_3_only_returns_one_of_two_sequences(self, tmpdir,
                                                                                estimated_count_threshold):
        # given
        kmer_size = 29
        seq_prefix = 'ATA' * 10

        b = FastqBuilder(tmpdir / 'single.fq')
        [b.with_seq(seq_prefix + 'ACAAC') for _ in range(10)]
        [b.with_seq(seq_prefix[2:] + 'ACAACCC') for _ in range(8)]
        [b.with_seq(seq_prefix[2:] + 'ACAACGG') for _ in range(4)]

        input_fastq = b.build()

        out_dir = Path(tmpdir) / 'abeona'
        args = [
            '--bootstrap-proportion-threshold', 1,
            '--estimated-count-threshold', estimated_count_threshold,
            '--fastx-single', input_fastq,
            '--kallisto-fragment-length', 4,
            '--kallisto-sd', 0.1,
            '--bootstrap-samples', 100,
            '--out-dir', out_dir,
            '--kmer-size', kmer_size,
            '--min-unitig-coverage', 0,
        ]

        # when
        AbeonaRunner().assemble(*args)

        # then
        expect = AbeonaExpectation(out_dir)
        if estimated_count_threshold == 5:
            expect.has_out_all_transcripts(seq_prefix + 'ACAACCC')
        elif estimated_count_threshold == 0:
            expect.has_out_all_transcripts(seq_prefix + 'ACAACCC', seq_prefix + 'ACAACGG')
        else:
            raise Exception


class TestCandidateTranscriptAnnotation:
    def test_estimates_abundances_of_two_transcripts_as_four(self, tmpdir):
        # given
        kmer_size = 29
        seq_prefix = 'ATA' * 10

        b = FastqBuilder(tmpdir / 'single.fq')
        [b.with_seq(seq_prefix + 'ACAAC') for _ in range(8)]
        [b.with_seq(seq_prefix[2:] + 'ACAACCC') for _ in range(4)]
        [b.with_seq(seq_prefix[2:] + 'ACAACGG') for _ in range(4)]

        input_fastq = b.build()

        out_dir = Path(tmpdir) / 'abeona'
        args = [
            '--estimated-count-threshold', 2,
            '--fastx-single', input_fastq,
            '--kallisto-fragment-length', 4,
            '--kallisto-sd', 0.1,
            '--bootstrap-samples', 100,
            '--out-dir', out_dir,
            '--kmer-size', kmer_size,
            '--min-unitig-coverage', 0,
        ]

        # when
        AbeonaRunner().assemble(*args)

        # then
        expect = AbeonaExpectation(out_dir)
        expect.has_out_all_transcripts(seq_prefix + 'ACAACCC', seq_prefix + 'ACAACGG')
        expect.has_out_all_transcript_description_matching(
            'prop_bs_est_counts_ge_2.0=[\d\.]+;est_count=8')


class TestBugsFromUsers:
    def test_ignores_subgraph_that_does_not_pseudoalign(self, tmpdir):
        # given
        b = PairedFastqBuilder(tmpdir / 'forward.fq', tmpdir / 'reverse.fq')
        read_pairs = list(
            zip(
                (l.lstrip() for l in """AGGACATAAAGCTCACGTGAGGTAATTACTTGAAGAAATTCCATTTAGAAGTTGTTCCTCTAGACATTGACTGGGCTGCCTAAGCCAACGTACTTCAATGGACTGAGGACATAAAGCTCACGTGAG
                GACATTGATTGAGCTGCCTAAGCCAACGTAGTTCAATGGATTGAGGACATAAAGCTCACGTGAGGTAATTACTTGAAGAAATTCCATTTAGAAGTTGTTCCTCTAGACATTGACTGGGCTGCCTAA
                GACATTGATTGAGCTGCCTAAGCCAACGTAGTTCAATGGATTGAGGACATAAAGCTCACGTGAGGTAATTACTTGAAGAAATTCCATTTAGAAGTTGTTCCTCTAGACATTGACTGGGCTGCCTAA
                GAGACATTGATTGAGCTGCCTAAGCCAACGTAGTTCAATGGATTGAGGACATAAAGCTCACGTGAGGTAATTACTTGAAGAAATTCCATTTAGAAGTTGTTCCTCTAGACATGGACTGGGCTGCCT
            """.split('\n')),
                (l.lstrip() for l in """CCGGGAGCAAATTGGGTTGAAGTCACAACAAAGTTTAAGCAATGGTATGAACAAGTTTTCGCGAGGATAGAGGAACAACTTCTAAATGGAATTTCTTCAAGTAATTACCTCACGTGAGCTTTATGT
                GCAAATTGGGTTGAAGTCACAACAAAGTTTAAGCAATGGTATGAACAAGTTTTCGCGAGGATAGAGGAACAACTTCTAAATGGAATTTCTTCAAGTAATTACCTCACGTGAGCTTTATGTCCTCAG
                GCAAATTGGGTTGAAGTCACAACAAAGTTTAAGCAATGGTATGAACAAGTTTTCGCGAGGATAGAGGAACAACTTCTAAATGGAATTTCTTCAAGTAATTACCTCACGTGAGCTTTATGTCCTCAG
                CACAACAAAGTTTAAGCAATGGTATGAACAAGTTTTCGCGAGGATAGAGGAACAACTTCTAAATGGAATTTCTTCAAGTAATTACCTCACGTGAGCTTTATGTCCTCAGTCCATTGAAGTACGTTG
            """.split('\n'))
            )
        )
        print(read_pairs)
        for rp in read_pairs:
            b.with_pair(*rp)

        forward, reverse = b.build()

        out_dir = Path(tmpdir) / 'abeona'
        args = [
            '--fastx-forward', forward,
            '--fastx-reverse', reverse,
            '--bootstrap-samples', 10,
            '--out-dir', out_dir,
            '--kmer-size', 55,
        ]

        # when
        AbeonaRunner().assemble(*args)

        # then
        expect = AbeonaExpectation(out_dir)
        expect.has_out_all_transcripts()

    def test_ignores_subgraphs_where_cortexpy_link_traversal_fails(self, tmpdir):
        # given
        fixtures = Path('tests/fixtures/2018-12-13')
        forward, reverse = fixtures / '1.fa', fixtures / '2.fa'
        out_dir = Path(tmpdir) / 'abeona'
        args = [
            '--fastx-forward', forward,
            '--fastx-reverse', reverse,
            '--bootstrap-samples', 10,
            '--out-dir', out_dir,
            '--kmer-size', 47,
        ]

        # when/then (no error)
        AbeonaRunner().assemble(*args)
