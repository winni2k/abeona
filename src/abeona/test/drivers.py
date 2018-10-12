import shutil
import subprocess
from pathlib import Path

import attr
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from cortexpy.test import builder

from abeona.cli import main as abeona_main
from .expectations import Traversals


@attr.s(slots=True)
class SubgraphTestDriver(object):
    tmpdir = attr.ib()
    builder = attr.ib(attr.Factory(builder.Mccortex))
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
class ReadsTestDriver(object):
    tmpdir = attr.ib()
    builder = attr.ib(attr.Factory(builder.Mccortex))

    def __getattr__(self, item):
        if item in ['with_kmer_size', 'with_dna_sequence']:
            return getattr(self.builder, item)
        else:
            raise ValueError(f'Could not find {item}')

    def run(self):
        graphs = SubgraphTestDriver(self.tmpdir, self.builder).build().traversals
        subprocess.run(f'cat {tmpdir}/*.fasta > {tmpdir}/combined.fasta', shell=True, check=True)

        with open(tmpdir / 'subgraph_list.txt', 'w') as fh:
            for graph in graphs:
                fh.write(f'{graph.with_suffix()}')

        command = f'abeona reads {graph} {tmpdir}/combined.fasta '
