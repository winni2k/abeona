import shutil
import subprocess
from collections import OrderedDict
from logging import getLogger
from pathlib import Path

import attr
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from cortexpy.test import runner

from abeona.cli import main as abeona_main
from .expectations import Traversals, Fastas, TraversalsAndFastas

logger = getLogger()


@attr.s(slots=True)
class Mccortex:
    kmer_size = attr.ib(3)
    sequences = attr.ib(attr.Factory(OrderedDict))
    sequence_pairs = attr.ib(attr.Factory(OrderedDict))
    mccortex_bin = attr.ib('mccortex')

    def assert_single_or_paired(self):
        if len(self.sequences) > 0 and len(self.sequence_pairs) > 0:
            raise NotImplementedError(
                'Specifying both single and paired end reads to driver is not implemented')

    @property
    def is_paired(self):
        self.assert_single_or_paired()
        assert len(self.sequence_pairs) != 0 or len(self.sequences) != 0
        return len(self.sequence_pairs) > 0

    def with_kmer_size(self, kmer_size):
        self.kmer_size = kmer_size
        return self

    def with_dna_sequence(self, sequence, *, name='sample_0'):
        if name not in self.sequences:
            self.sequences[name] = []
        self.sequences[name].append([sequence])
        return self

    def with_dna_sequence_pair(self, sequence1, sequence2, name='sample_0'):
        if name not in self.sequence_pairs:
            self.sequence_pairs[name] = []
        self.sequence_pairs[name].append([sequence1, sequence2])
        return self

    def build(self, dir):
        dir = Path(dir)
        dir.mkdir(exist_ok=True)
        self.assert_single_or_paired()
        if len(self.sequence_pairs) == 0:
            sequences = self.sequences
        else:
            sequences = self.sequence_pairs

        mccortex_args = f'build --force --sort --kmer {self.kmer_size}'
        for name, dna_sequence_tuple_list in sequences.items():
            mccortex_args += f' --sample {name}'
            is_paired = len(dna_sequence_tuple_list[0]) == 2
            fa1 = dir / f'input.{name}.1.fasta'
            fa2 = dir / f'input.{name}.2.fasta'
            with open(fa1, 'w') as fh1:
                with open(fa2, 'w') as fh2:
                    for tuple in dna_sequence_tuple_list:
                        fh1.write(SeqRecord(Seq(tuple[0]), id='0', description='').format('fasta'))
                        if is_paired:
                            fh2.write(
                                SeqRecord(Seq(tuple[1]), id='0', description='').format('fasta'))
            if is_paired:
                mccortex_args += f' -2 {fa1}:{fa2}'
            else:
                mccortex_args += f' -1 {fa1}'

        output_graph = str(dir / 'output.ctx')
        mccortex_args += f' {output_graph}'

        ret = runner.Mccortex(self.kmer_size, mccortex_bin=self.mccortex_bin) \
            .run(mccortex_args.split())
        logger.debug('\n' + ret.stdout.decode())
        logger.debug('\n' + ret.stderr.decode())

        ret = runner.Mccortex(self.kmer_size, mccortex_bin=self.mccortex_bin) \
            .view(output_graph)
        logger.debug('\n' + ret.stdout.decode())

        return output_graph


@attr.s(slots=True)
class CleanMccortex:
    mc_builder = attr.ib()
    min_unitig_coverage = attr.ib(0)
    min_tip_length = attr.ib(0)

    @classmethod
    def from_mccortex_builder(cls, builder):
        return cls(mc_builder=builder)

    def __getattr__(self, item):
        return getattr(self.mc_builder, item)

    def with_min_unitig_coverage(self, n):
        self.min_unitig_coverage = n
        return self

    def build(self, dir):
        dir = Path(dir)
        dir.mkdir(exist_ok=True)
        full_graph = self.mc_builder.build(dir)
        clean_graph = dir / 'clean.ctx'

        ret = runner.Mccortex(self.mc_builder.kmer_size, mccortex_bin=self.mc_builder.mccortex_bin) \
            .clean(graphs=[full_graph], out=clean_graph, clip_tips_shorter_than=self.min_tip_length,
                   unitigs_with_mean_coverage_less_than=self.min_unitig_coverage)
        logger.debug('\n' + ret.stdout.decode())
        return clean_graph

@attr.s(slots=True)
class SubgraphTestDriver:
    tmpdir = attr.ib()
    builder = attr.ib(attr.Factory(Mccortex))
    initial_contigs = attr.ib(None)

    def __getattr__(self, item):
        if item in ['with_kmer_size', 'with_dna_sequence']:
            return getattr(self.builder, item)
        else:
            raise ValueError(f'Could not find {item}')

    def with_initial_contigs(self, *contigs):
        self.initial_contigs = contigs
        return self

    def run(self):
        shutil.rmtree(self.tmpdir)
        self.tmpdir.mkdir()
        graph = self.builder.build(self.tmpdir)
        out_dir = self.tmpdir / 'abeona_output'
        command = ['subgraphs', graph, out_dir]

        if self.initial_contigs:
            initial_contigs = self.tmpdir / 'initial-kmers.fa'
            records = [
                SeqRecord(Seq(contig), id=f'initial_contig_{idx}') for idx, contig in
                enumerate(self.initial_contigs)
            ]
            SeqIO.write(records, str(initial_contigs), 'fasta')
            command += ['--initial-contigs', initial_contigs]

        args = ['abeona'] + [str(arg) for arg in command]
        abeona_main(args)

        graphs = list(Path(out_dir).glob('g*.traverse.ctx'))

        return Traversals(graphs)


@attr.s(slots=True)
class ReadsTestDriver:
    builder = attr.ib(attr.Factory(Mccortex))
    record_buffer_size = attr.ib(None)

    def __getattr__(self, item):
        if item in ['with_kmer_size', 'with_dna_sequence', 'with_dna_sequence_pair', 'with_min_unitig_coverage']:
            return getattr(self.builder, item)
        else:
            raise ValueError(f'Could not find {item}')

    def with_record_buffer_size(self, n):
        self.record_buffer_size = n
        return self

    def run(self, tmpdir):
        out_dir = tmpdir / 'abeona_reads'
        traversal_expectation = SubgraphTestDriver(tmpdir, self.builder).run()
        graphs = traversal_expectation.traversals

        for idx in [1, 2]:
            subprocess.run(
                f'cat {tmpdir}/*.{idx}.fasta > {tmpdir}/combined.{idx}.fasta',
                shell=True,
                check=True
            )

        graph_list = tmpdir / 'subgraph_list.txt'
        with open(graph_list, 'w') as fh:
            fh.write('# prefix\tgraph\n')
            for graph in graphs:
                fh.write(f'{out_dir/graph.stem}\t{graph}\n')

        command = f'abeona reads --format fasta {graph_list} {tmpdir}/combined.1.fasta'
        if self.builder.is_paired:
            command += f' --reverse {tmpdir}/combined.2.fasta'
        if self.record_buffer_size:
            command += f' --record-buffer-size {self.record_buffer_size}'
        abeona_main(command.split())

        reads = list(Path(out_dir).glob('g*.fa.gz'))

        return TraversalsAndFastas(Fastas(reads), traversal_expectation)
