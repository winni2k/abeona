import attr
import networkx as nx
from cortexpy.graph.parser import RandomAccess
from cortexpy.graph.parser.streaming import kmer_list_generator_from_stream, load_cortex_graph
from cortexpy.graph.traversal import Engine
from cortexpy.test.expectation import KmerGraphExpectation
from cortexpy.utils import lexlo


@attr.s(slots=True)
class Graph(object):
    mccortex_graph = attr.ib()
    traversal = attr.ib(init=False)
    kmers = attr.ib(init=False)
    kmer_strings = attr.ib(init=False)

    def __attrs_post_init__(self):
        with open(self.mccortex_graph, 'rb') as fh:
            self.kmer_strings = set([''.join(k) for k in kmer_list_generator_from_stream(fh)])
            self.traversal = Engine(ra_parser=RandomAccess(fh)) \
                .traverse_from_each_kmer_in_iterable(self.kmer_strings) \
                .graph

    def has_n_subgraphs(self, n):
        subgraphs = list(nx.weakly_connected_component_subgraphs(self.traversal))
        assert n == len(subgraphs)
        return self

    def has_kmer_strings(self, *kmer_strings):
        assert self.kmer_strings == set(kmer_strings)
        return self


@attr.s(slots=True)
class Graphs(object):
    mccortex_graphs = attr.ib()
    graph_expectations = attr.ib(init=False)

    def __attrs_post_init__(self):
        self.graph_expectations = [Graph(graph) for graph in self.mccortex_graphs]

    def has_graph_with_kmers(self, *kmers):
        kmers = set(kmers)
        for expect in self.graph_expectations:
            if kmers == expect.kmer_strings:
                return self
        assert False

    def has_n_graphs(self, n):
        assert n == len(self.graph_expectations)
        return self


@attr.s(slots=True)
class Traversals(object):
    traversals = attr.ib()
    traversal_expectations = attr.ib(init=False)

    def __attrs_post_init__(self):
        self.traversal_expectations = [
            KmerGraphExpectation(load_cortex_graph(open(graph, 'rb'))) for graph in
            self.traversals
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
