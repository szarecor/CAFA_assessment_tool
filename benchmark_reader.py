import sys
from collections import defaultdict
import numpy
import os
from Ontology.IO import OboIO

legal_species = [
    "all",
    "eukarya",
    "prokarya",
    'HELPY',
    'ECOLI',
    'RAT',
    'DANRE',
    'SULSO',
    'DROME',
    'PSEPK',
    'STRPN',
    'PSEAE',
    'BACSU',
    'MYCGE',
    'HUMAN',
    'METJA',
    'DICDI',
    'YEAST',
    'SCHPO',
    'ARATH',
    'XENLA',
    'MOUSE',
    'PSESM',
    'SALTY',
    'CANAX',
    'SALCH']
legal_types = ["type1","type2","all"]
legal_subtypes = ["easy","hard"]
#easy and hard are only in NK benchmarks!!!!
legal_modes = ["full", "partial"]
root_terms = ['GO:0008150','GO:0005575','GO:0003674']


def go_ontology_ancestors_split_write(obo_path):
    """
    Input: an OBO file
    Output: 3 files with ancestors
    by Dr. Friedberg
    Updated 20190906 by Ashley: ony write ancestor within the same ontology
    """

    obo_bpo_out = open("%s_ancestors_bpo.txt" % (os.path.splitext(obo_path)[0]), "w")
    obo_cco_out = open("%s_ancestors_cco.txt" % (os.path.splitext(obo_path)[0]), "w")
    obo_mfo_out = open("%s_ancestors_mfo.txt" % (os.path.splitext(obo_path)[0]), "w")
    obo_parser = OboIO.OboReader(open(obo_path))
    go = obo_parser.read()
    mfo_terms, bpo_terms, cco_terms = go_ontology_split(go)
    for term in mfo_terms:
        ancestors = go.get_ancestors(term)
        ancestors_filter = set()
        if len(ancestors) > 0:
            for anc in ancestors:
                anc_ont = go.get_namespace(anc)
                if anc_ont == 'molecular_function':
                    ancestors_filter.add(anc)
        obo_mfo_out.write("%s\t%s\n" % (term, ",".join(ancestors_filter)))
    for term in bpo_terms:
        ancestors = go.get_ancestors(term)
        ancestors_filter = set()
        if len(ancestors) > 0:
            for anc in ancestors:
                anc_ont = go.get_namespace(anc)
                if anc_ont == 'biological_process':
                    ancestors_filter.add(anc)
        obo_bpo_out.write("%s\t%s\n" % (term, ",".join(ancestors_filter)))
    for term in cco_terms:
        ancestors = go.get_ancestors(term)
        ancestors_filter = set()
        if len(ancestors) > 0:
            for anc in ancestors:
                anc_ont = go.get_namespace(anc)
                if anc_ont == 'cellular_component':
                    ancestors_filter.add(anc)
        obo_cco_out.write("%s\t%s\n" % (term, ",".join(ancestors_filter)))

    obo_mfo_out.close()
    obo_bpo_out.close()
    obo_cco_out.close()
    return ([len(bpo_terms), len(cco_terms), len(mfo_terms)])


def go_ontology_split(ontology):
    """
    Split an GO obo file into three ontologies
    by Dr. Friedberg
    """
    mfo_terms = set({})
    bpo_terms = set({})
    cco_terms = set({})

    for node in ontology.get_ids():
        # loop over node IDs and alt_id's
        if ontology.namespace[node] == "molecular_function":
            mfo_terms.add(node)
        elif ontology.namespace[node] == "biological_process":
            bpo_terms.add(node)
        elif ontology.namespace[node] == "cellular_component":
            cco_terms.add(node)
        else:
            raise(ValueError, "{} has no namespace".format(node))

    return mfo_terms, bpo_terms, cco_terms


def read_benchmark(namespace, species, types, fullbenchmarkfolder, obopath):
    '''
    Read Benchmark.

    Input:
    namespace
    species
    types
    fullbenchmarkfolder
    obopath

    Output:
    bench
    obocountDict
    '''
    # Ancestor files here are precomputed
    # To get the ancestor files, use preprocess.py (go_ontology_ancestors_split_write) DOES NOT EXIST
    legal_types = ["type1", "type2", "typex"]
    legal_subtypes = ["easy", "hard"]
    legal_namespace = ["bpo", "mfo", "cco", "hpo"]
    # fullbenchmarkfolder = './precrec/benchmark/'
    if namespace not in legal_namespace:
        sys.stderr.write("Namespace not accepted, choose from 'bpo', 'cco', 'mfo' and 'hpo'\n")
    elif (species not in legal_species) and (species not in legal_subtypes):
        sys.stderr.write('Species not accepted')
    elif types not in legal_types:
        sys.stderr.write('Type not accepted, choose from "type1","type2" and "typex"\n')
    else:
        matchname = namespace + '_' + species + '_' + types + '.txt'
    # generate ancestor files
    obocounts = go_ontology_ancestors_split_write(obopath)
    obocountDict = {'bpo': obocounts[0], 'cco': obocounts[1], 'mfo': obocounts[2]}
    # ontology-specific calculations
    if namespace == 'bpo':
        full_benchmark_path = fullbenchmarkfolder + '/groundtruth/' + 'leafonly_BPO.txt'
        ancestor_path = os.path.splitext(obopath)[0] + "_ancestors_bpo.txt"
    elif namespace == 'cco':
        full_benchmark_path = fullbenchmarkfolder + '/groundtruth/' + 'leafonly_CCO.txt'
        ancestor_path = os.path.splitext(obopath)[0] + "_ancestors_cco.txt"
    elif namespace == 'mfo':
        full_benchmark_path = fullbenchmarkfolder + '/groundtruth/' + 'leafonly_MFO.txt'
        ancestor_path = os.path.splitext(obopath)[0] + "_ancestors_mfo.txt"
    benchmarkListPath = fullbenchmarkfolder + '/lists/' + matchname
    if os.path.isfile(benchmarkListPath) and os.path.getsize(benchmarkListPath) > 0:
        handle = open(fullbenchmarkfolder + '/lists/' + matchname, 'r')
        prots = set()
        for line in handle:
            prots.add(line.strip())
        handle.close()
        tempfilename = 'temp_%s_%s_%s.txt' % (namespace, species, types)
        tempfile = open(fullbenchmarkfolder + '/' + tempfilename, 'w')
        for line in open(full_benchmark_path, 'r'):
            prot = line.split('\t')[0]
            if prot in prots:
                tempfile.write(line)
        tempfile.close()
        bench = benchmark(ancestor_path, tempfile.name)
        bench.propagate()
        os.remove(tempfile.name)
    else:
        print('Benchmark set is empty.\n')
        bench = None
    return bench, obocountDict


class benchmark:
    def __init__(self, ancestor_path, benchmark_path):
        '''
        Initialize the benchmark.

        Input:
        benchmark_path is ontology specific
        ancestor_path is ontology specific
        '''

        # Key: protein
        # Value: set of benchmark leaf terms
        self.ancestors = defaultdict(set)
        # Read GO ancestors file generated with go_ontology_ancestors_split_write()
        # File format:
        # go_term <tab> ancestor_1,ancestor_2,..,ancestor_n
        with open(ancestor_path) as ancestors_input:
            for inline in ancestors_input:
                inrec = inline.strip().split('\t')
                term = inrec[0]
                if len(inrec) == 1:
                    self.ancestors[term] = set({})
                else:
                    term_ancestors = inrec[1]
                    self.ancestors[term] = set(term_ancestors.split(','))

        self.true_base_terms = defaultdict(set)
        with open(benchmark_path) as benchmark_input:
            for inline in benchmark_input:
                protein, term = inline.strip().split('\t')
                self.true_base_terms[protein].add(term)

    def propagate(self):
        '''
        Progate Benchmark terms.
        '''

        # Key: protein
        # Value: set of benchmark propagated terms
        self.true_terms = defaultdict(set)
        for protein in self.true_base_terms:
            for term in self.true_base_terms[protein]:
                try:

                    ancestors = self.ancestors[term].difference(root_terms)
                    # delete root term in self.true_terms
                # modified on 20170203
                except KeyError:
                    sys.stderr.write("%s not found \n" % term)
                self.true_terms[protein].add(term)
                self.true_terms[protein] |= ancestors


if __name__ == "__main__":
    '''
   bpo HUMAN type1 
 
    legal_types = ["type1", "type2", "typex"]
    legal_subtypes = ["easy", "hard"]
    legal_namespace = ["bpo", "mfo", "cco", "hpo"]
    '''

    namespace = "bpo"
    species = "HUMAN"
    types = "type1"
    benchmark_folder = "./data/benchmark/CAFA3_benchmarks/"
    obopath = "./data/go_cafa3.obo"
    legal_types = ["type1", "type2", "typex"]
    legal_subtypes = ["easy", "hard"]
    legal_namespace = ["bpo", "mfo", "cco", "hpo"]

    benchmark, obo_counts = read_benchmark(namespace, species, types, benchmark_folder, obopath)
    print("BENCHMARK:")
    print(benchmark)
    print("---------------------------------")
    print("OBO COUNTS:", obo_counts)
