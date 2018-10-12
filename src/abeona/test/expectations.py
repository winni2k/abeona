import gzip
from pathlib import Path

import attr
from cortexpy.graph.parser.streaming import load_cortex_graph
from cortexpy.test.expectation import KmerGraphExpectation, Fasta
from cortexpy.utils import lexlo


# @attr.s(slots=True)
# class Graph:
#     mccortex_graph = attr.ib()
#     traversal = attr.ib(init=False)
#     kmers = attr.ib(init=False)
#     kmer_strings = attr.ib(init=False)
#
#     def __attrs_post_init__(self):
#         with open(self.mccortex_graph, 'rb') as fh:
#             self.kmer_strings = set([''.join(k) for k in kmer_list_generator_from_stream(fh)])
#             self.traversal = Engine(ra_parser=cortexpy.graph.parser.RandomAccess(fh)) \
#                 .traverse_from_each_kmer_in_iterable(self.kmer_strings) \
#                 .graph
#
#     def has_n_subgraphs(self, n):
#         subgraphs = list(nx.weakly_connected_component_subgraphs(self.traversal))
#         assert n == len(subgraphs)
#         return self
#
#     def has_kmer_strings(self, *kmer_strings):
#         assert self.kmer_strings == set(kmer_strings)
#         return self


# @attr.s(slots=True)
# class Graphs:
#     mccortex_graphs = attr.ib()
#     graph_expectations = attr.ib(init=False)
#
#     def __attrs_post_init__(self):
#         self.graph_expectations = [Graph(graph) for graph in self.mccortex_graphs]
#
#     def has_graph_with_kmers(self, *kmers):
#         kmers = set(kmers)
#         for expect in self.graph_expectations:
#             if kmers == expect.kmer_strings:
#                 return self
#         assert False
#
#     def has_n_graphs(self, n):
#         assert n == len(self.graph_expectations)
#         return self
#

@attr.s(slots=True)
class AbeonaKmerGraphExpectation(KmerGraphExpectation):
    graph_path = attr.ib(None)


@attr.s(slots=True)
class Traversals:
    traversals = attr.ib()
    traversal_expectations = attr.ib(init=False)

    def __attrs_post_init__(self):
        self.traversal_expectations = [
            AbeonaKmerGraphExpectation(load_cortex_graph(open(graph, 'rb')), graph_path=graph) for
            graph in self.traversals
        ]

    def has_graph_with_kmers(self, *kmers):
        kmers = set(kmers)
        for expect in self.traversal_expectations:
            if kmers == set(lexlo(kmer_string) for kmer_string in expect.graph):
                return expect
        assert False

    def has_n_graphs(self, n):
        assert n == len(self.traversal_expectations)
        return self


@attr.s(slots=True)
class Fastas:
    fastas = attr.ib()

    def has_fasta_with_prefix(self, stem):
        for fasta in self.fastas:
            if Path(fasta).with_suffix('').stem == stem:
                with gzip.open(fasta, 'rt') as fh:
                    return Fasta(fh.read())
        assert False, f'Could not find stem ({stem}) in fastas\n{self.fastas}'


@attr.s(slots=True)
class TraversalsAndFastas:
    fastas_expection = attr.ib()
    traversals_expectation = attr.ib()

    def has_graph_with_kmers_and_mapped_reads(self, kmers, seqs):
        graph = self.traversals_expectation.has_graph_with_kmers(*kmers).graph_path
        graph_stem = Path(graph).stem

        fasta_expect = self.fastas_expection.has_fasta_with_prefix(graph_stem)
        for seq in seqs:
            fasta_expect.has_record(seq)
        fasta_expect.has_n_records(len(seqs))
        return self

    def has_n_graphs(self,n):
        self.traversals_expectation.has_n_graphs(n)
        return self
