import json
from pathlib import Path
from contextlib import closing

import attr
from cortexpy.graph.parser.streaming import load_cortex_graph
from cortexpy.test.expectation import KmerGraphExpectation, Fasta
from .utils import get_maybe_gzipped_file_handle

@attr.s(slots=True)
class AbeonaKmerGraphExpectation(KmerGraphExpectation):
    graph_path = attr.ib(None)
    meta = attr.ib(None)

    def has_meta_info(self, key, val):
        assert val == self.meta[key]
        return self


@attr.s(slots=True)
class Traversals:
    traversals = attr.ib()
    traversal_expectations = attr.ib(init=False)

    def __attrs_post_init__(self):
        self.traversal_expectations = [
            AbeonaKmerGraphExpectation(
                load_cortex_graph(open(graph, 'rb')),
                graph_path=graph,
                meta=json.load(open(str(graph) + '.json', 'r'))
            ) for graph in self.traversals
        ]

    def has_graph_with_kmers(self, *kmers):
        kmers = set(kmers)
        tested_kmers = []
        for expect in self.traversal_expectations:
            tested_kmers.append(sorted(expect.graph))
            if sorted(kmers) == tested_kmers[-1]:
                return expect
        assert False, f'Could not find graph with kmers: {kmers}\nin: {tested_kmers}'

    def has_n_graphs(self, n):
        assert n == len(self.traversal_expectations)
        return self


@attr.s(slots=True)
class Fastas:
    fastas = attr.ib()

    def has_fastas_with_prefix(self, stem):
        fastas = []
        for fasta in sorted(self.fastas):
            if Path(fasta).with_suffix('').stem.startswith(stem):
                with closing(get_maybe_gzipped_file_handle(fasta, 'rt')) as fh:
                    fastas.append(Fasta(fh.read()))
        if len(fastas) == 0:
            assert False, f'Could not find stem ({stem}) in fastas\n{self.fastas}'
        return fastas


@attr.s(slots=True)
class TraversalsAndFastas:
    fastas_expection = attr.ib()
    traversals_expectation = attr.ib()

    def has_graph_with_kmers_and_mapped_reads(self, kmers, seqs, seqs2=None):
        graph = self.traversals_expectation.has_graph_with_kmers(*kmers).graph_path
        graph_stem = Path(graph).stem

        fasta_expects = self.fastas_expection.has_fastas_with_prefix(graph_stem)

        seq_lists = [seqs]
        if seqs2 is not None:
            seq_lists.append(seqs2)
        assert len(fasta_expects) == len(seq_lists), f'len({fasta_expects}) != len({seq_lists})'
        for idx, seq_list in enumerate(seq_lists):
            for seq in seq_list:
                fasta_expects[idx].has_record(seq)
            fasta_expects[idx].has_n_records(len(seq_list))
        return self

    def has_n_graphs(self, n):
        self.traversals_expectation.has_n_graphs(n)
        return self
