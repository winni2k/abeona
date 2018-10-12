import shutil
import subprocess
from pathlib import Path

import attr
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from cortexpy.test import builder

from abeona.cli import main as abeona_main
from .expectations import Traversals, Fastas, TraversalsAndFastas


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
        out_dir = self.tmpdir / 'abeona_reads'
        traversal_expectation = SubgraphTestDriver(self.tmpdir, self.builder).run()
        graphs = traversal_expectation.traversals

        subprocess.run(f'cat {self.tmpdir}/*.fasta > {self.tmpdir}/combined.fasta', shell=True, check=True)

        graph_list = self.tmpdir / 'subgraph_list.txt'
        with open(graph_list, 'w') as fh:
            fh.write('# prefix\tgraph\n')
            for graph in graphs:
                fh.write(f'{out_dir/graph.stem}\t{graph}\n')

        command = f'abeona reads --format fasta {graph_list} {self.tmpdir}/combined.fasta'
        abeona_main(command.split())

        reads = list(Path(out_dir).glob('g*.fasta.gz'))

        return TraversalsAndFastas(Fastas(reads), traversal_expectation)
