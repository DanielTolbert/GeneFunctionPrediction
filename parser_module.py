import networkx as nx
import re
import pandas as pd
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
from operator import add

import warnings
warnings.filterwarnings("ignore")

NAMESPACES = ['biological_process', 'molecular_function', 'cellular_component']

ROOT_TERMS = dict.fromkeys(NAMESPACES, None)

ASPECTS = {'biological_process': 'P', 'molecular_function': 'F', 'cellular_component': 'C'}


CODE_GROUPS = {'EXP': {'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP'}, 'HTP': {'HTP', 'HDA', 'HMP', 'HGI', 'HEP'}, 'PHYLO': {'IBA', 'IBD', 'IKR', 'IRD'}, 'COMP': {'ISS', 'ISO', 'ISA', 'ISM', 'IGC', 'RCA'}, 'AUTH': {'TAS', 'NAS', 'IC'}, 'IEA': {'IEA'}, 'ND': {'ND'}}

class Annotation_Parser:
    def __init__(self, obo_file, gaf_file, threshold=1):
        self._threshold = threshold
        if obo_file is not None:
            self._init_GO_graph(obo_file)
        if gaf_file is not None:
            self._init_annotations(gaf_file)
        self._S = self.calculate_S()
        self._T = self.calculate_T()

    @property
    def GO_graph(self) -> nx.DiGraph():
        return self._GO_graph

    @GO_graph.setter
    def GO_graph(self, new_graph: nx.DiGraph):
        self._GO_graph = new_graph

    @property
    def threshold(self):
        return self._threshold

    @threshold.setter
    def threshold(self, new_threshold):
        self._threshold = new_threshold
        self.S = self.calculate_S()
        self.T = self.calculate_T()
        
    @property
    def ancestors_by_GO(self):
        return self._ancestors_by_GO

    @ancestors_by_GO.setter
    def ancestors_by_GO(self, new_ancestors):
        self._ancestors_by_GO = new_ancestors
    
    @property
    def annotations(self):
        return self._annotations

    @annotations.setter
    def annotations(self, new_annotations):
        self._annotations = new_annotations

    """
    all GO terms that annotate at least a threshold percent of genes
    """
    @property
    def S(self):
        return self._S
    
    @S.setter
    def S(self, new_S):
        self._S = new_S

    """
    most specific GO terms that annotate at least a threshold percent of genes
    """
    @property
    def T(self):
        return self._T
    
    @T.setter
    def T(self, new_T):
        self._T = new_T

    
    @property
    def n_bp(self):
        return self._n_bp

    @n_bp.setter
    def n_bp(self, new_n):
        self._n_bp = new_n

    @property
    def n_mf(self):
        return self._n_mf

    @n_mf.setter
    def n_mf(self, new_n):
        self._n_mf = new_n

    @property
    def n_cc(self):
        return self._n_cc

    @n_cc.setter
    def n_cc(self, new_n):
        self._n_cc = new_n
    
    def _init_GO_graph(self, obo_file):
        G = nx.DiGraph()
        with open(obo_file, 'r') as input:
            for line in input:
                if line.startswith('[Term]'):
                    id = next(input, '').strip('\n').split(':', 1)[1].strip()
                    name = next(input, '').strip('\n').split(':', 1)[1].strip()
                    namespace = next(input, '').strip('\n').split(':', 1)[1].strip()
                    line2 = next(input, '').strip('\n')
                    G.add_node(id, name=name, namespace=namespace)
                    while line2 != '':
                        if line2.startswith('is_obsolete'):
                            IS_OBSOLETE = True
                            G.remove_node(id)
                            break
                        elif line2.startswith('is_a'):
                            IS_A_ENCOUNTERED = True
                            contents = line2.split(':', 1)[1]
                            other = contents.split('!')[0].strip()
                            G.add_edge(id, other, type='is_a')
                        elif line2.startswith('relationship'):
                            contents = line2.split(':', 1)[1].split('!', 1)[0]
                            rel_type_other = contents.split('GO')
                            rel_type = rel_type_other[0].strip()
                            other = 'GO'+rel_type_other[1].strip()
                            G.add_edge(id, other, type=rel_type)
                        line2 = next(input, '').strip('\n')
                    
                    # if not IS_A_ENCOUNTERED and not IS_OBSOLETE:
                    #     print(id)

                    IS_OBSOLETE = False
                    IS_A_ENCOUNTERED = False
        self.GO_graph = G

    def _init_annotations(self, gaf_file):
        df = pd.read_csv(gaf_file, comment='!', sep='\t', header=None, names=['DB', 'DB_Object_ID', 'DB_Object_Symbol', 'Qualifier', 'GO_ID', 'DB:Reference', 'Evidence_Code', 'With_or_From', 'Aspect', 'DB_Object_Name', 'DB_Object_Synonym', 'DB_Object_Type', 'Taxon', 'Date', 'Assigned_By', 'Annotation_Extension', 'Gene_Product_Form_ID'])

        df = df[~df.Qualifier.str.contains("NOT")]

        d = df[['DB_Object_ID', 'Evidence_Code', 'Aspect']].agg(tuple, 1).groupby(df['GO_ID']).apply(list).to_dict()

        ancestors = {}
        human_annotations_transferred = defaultdict(list)
        for i, u in enumerate(self.GO_graph.nodes()):
            ancestors[u] = nx.descendants(self.GO_graph, u)
            if len(ancestors[u]) == 0:
                name = self.GO_graph.nodes[u]['name']
                ROOT_TERMS[name] = u
            try:
                annotated_genes = d[u]
            except Exception:
                continue

            for a in ancestors[u]:
                for annotation in annotated_genes:
                    gene = annotation[0]
                    evidence_code = annotation[1]
                    aspect = annotation[2]
                    human_annotations_transferred[a].append((gene, evidence_code, aspect))

        for namespace, GO in ROOT_TERMS.items():
            terms_annotated = len(human_annotations_transferred[GO])
            #print(namespace, terms_annotated)
            if namespace == 'biological_process':
                self.n_bp = terms_annotated
            elif namespace == 'molecular_function':
                self.n_mf = terms_annotated
            elif namespace == 'cellular_component':
                self.n_cc = terms_annotated

        self.annotations = human_annotations_transferred
        self.ancestors_by_GO = ancestors

    def calculate_S(self):
        THRESHOLD_PERCENT = self.threshold
        S = defaultdict(list)
        for GO in self.annotations:
            namespace = self.GO_graph.nodes[GO]['namespace']
            terms_annotated_by_GO = len(self.annotations[GO])

            if namespace == 'biological_process':
                if (terms_annotated_by_GO >= (self.n_bp*THRESHOLD_PERCENT/100)):
                    S[namespace].append(GO)
            elif namespace == 'molecular_function':
                if (terms_annotated_by_GO >= (self.n_mf*THRESHOLD_PERCENT/100)):
                    S[namespace].append(GO)
            elif namespace == 'cellular_component':
                if (terms_annotated_by_GO >= (self.n_cc*THRESHOLD_PERCENT/100)):
                    S[namespace].append(GO)
        return S
       

    def calculate_T(self):
        all_ancestors = defaultdict(list)
        T = defaultdict(list)
        for namespace, GO in ROOT_TERMS.items():
            for GO in self.S[namespace]:
                all_ancestors[namespace] += self.ancestors_by_GO[GO]

            for GO in self.S[namespace]:
                if GO not in all_ancestors[namespace]:
                    T[namespace].append(GO)
        return T

    def annotations_at_threshold(self):
        annotated_genes = defaultdict(dict)
        for namespace, GO in ROOT_TERMS.items():
            annotated_genes[namespace] = {}
            namespace_abbrev = ASPECTS[namespace]
            for i, GO in enumerate(self.T[namespace]):
                #print(len(self.annotations[GO]))
                for annotation in self.annotations[GO]:
                    aspect = annotation[2]
                    if aspect != namespace_abbrev:
                        continue
                    gene = annotation[0]
                    evidence_code = annotation[1]
                    if gene not in annotated_genes[namespace]:

                        annotated_genes[namespace][gene] = {evidence_code}
                    else:
                        annotated_genes[namespace][gene].add(evidence_code)
        return annotated_genes

    def annotations_by_evidence(self):
        annotated_genes = self.annotations_at_threshold()
        annotations_by_evidence = defaultdict(dict)
        for namespace, GO in ROOT_TERMS.items():
            annotations_by_evidence[namespace] = defaultdict(list, {k:[] for k in CODE_GROUPS.keys()})
            for gene in annotated_genes[namespace]:
                gene_codes = annotated_genes[namespace][gene]
                for group, group_codes in CODE_GROUPS.items():
                    if not group_codes.isdisjoint(gene_codes):
                        annotations_by_evidence[namespace][group].append(gene)
                        break
        return annotations_by_evidence

    def stacked_plot(self):
        annotations_by_evidence = self.annotations_by_evidence()

        counts = defaultdict(list, {k:[] for k in CODE_GROUPS.keys()})
        for namespace, GO in ROOT_TERMS.items():
            for group in annotations_by_evidence[namespace]:
                count = len(annotations_by_evidence[namespace][group])
                counts[group].append(count)

        for i, (group, count) in enumerate(counts.items()):
            if i == 0:
                plt.bar(NAMESPACES, count, label=group)
                prev_count = count
            else:
                plt.bar(NAMESPACES, count, bottom=prev_count, label=group)
                prev_count = list(map(add, prev_count, count) )

        MAX_COUNT = np.max(prev_count)

        plt.xlabel('Aspect')
        plt.ylabel('Number of genes annotated')
        plt.title('Number of annotations by evidence')
        plt.legend()
        plt.ylim((0,MAX_COUNT*1.1))
        plt.show()
