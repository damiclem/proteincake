import json
import gzip
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import fisher_exact
from wordcloud import WordCloud

"""
Perform Fisher test. An Odd-Ratio above 77 tells us the GO prefers the first dataframe (p-value < 0.05),
while an Odd-Ratio under 0.013 tells us the GO prefers the second dataframe.
Return a Dataframe with index the GO and values the Odd-Ratio and the p-value.
"""
def fisher_test(df1, df2, col_name_go = 'go_id'):

        # Inint dict
        results = {}

        # Get the number of occurrances of the GO counts
        dict1, dict2 = dict(df1[col_name_go].value_counts()), dict(df2[col_name_go].value_counts())

        # Compute the intersaction of the GO terms
        key_intersection = set(dict1.keys()).intersection(set(dict2.keys()))

        for key in key_intersection:
        ### 1. Set frequencies
        # Number of occurrences of the specific GO term in DF1
            tp = dict1[key]
            # Number of occurrences of the specific GO term in DF2
            tn = dict2[key]
            # Number of GO terms that are different from the specific one in DF1
            fp = sum(dict1.values()) - tp
            # Number of GO terms that are different from the specific one in DF2
            fn = sum(dict2.values()) - tn
            # 2. Perform Fisher Exact Test
            fisher_results = fisher_exact([[tp, tn],[fp, fn]])
            # 3. Save results
            results.setdefault(key, {'OddRatio': fisher_results[0], 'p-value': fisher_results[1]})

        # Return the DataFrame
        return pd.DataFrame(results).transpose()

"""
To retrieve the GO that are parents, we cycle over ontology["graphs"][0]["edges"] which is a list of dictionary.
Every dictionary is about a GO id (sub) with a relation (pred) with another GO (obj).
We create a dictionary (parents) with as keys the sons and as values the parents.
"""
def get_parents(ontology):
    parents = {}  # {GO_id(son) : list of GO_id (parents)}
    for edge in ontology["graphs"][0]["edges"]:
        # select only is_a edges
        if edge["pred"] == "is_a":
            parents.setdefault(edge["sub"].split("_")[1], []).append(edge["obj"].split("_")[1])
    return parents

"""
Here we cycle over the nodes to obtain a dictionary of GO_id with as value a description.
- ontology["graphs"][0]["nodes"] is a list of dictionary with dict_keys(['id', 'meta', 'type', 'lbl'])
- ontology["graphs"][0]["nodes"][1]['lbl'] is the value (e.g: "endocytosed synaptic vesicle processing via endosome")
"""
def get_labels(ontology):
    labels = {}  # {term (GO_id): definition}
    for node in ontology["graphs"][0]["nodes"]:
        # exclude obsolete terms
        if "GO_" in node["id"] and "deprecated" not in node["meta"]:
            go_id = node["id"].split("_")[1]
            labels[go_id] = node["lbl"]
    return labels

"""
Build an ancestors dictionary with as key an GO_id and as value a list of GO_id which are the ancestors of the key.
Return ancestors = {GO_id : list of ancestor GO_ids}
"""
def get_ancestors(ontology):
    nodes = get_labels(ontology).keys()
    parents = get_parents(ontology)
    ancestors = {}
    for node in nodes:
        node_ancestors = []
        node_parents = parents.get(node)
        # Loop parent levels until no more parents
        while node_parents:
            node_ancestors.extend(node_parents)
            # Get the parents of current parents (1 level up)
            node_parents = [term for parent in node_parents for term in parents.get(parent, [])]
        ancestors[node] = node_ancestors
    return ancestors

"""
Build a dictionary for the children (similar to the ancestors one)
Return {node : list_of_children}, leaf terms are not keys.
"""
def get_children(ontology):
    ancestors = get_ancestors(ontology)
    children = {}
    for node in ancestors:
        for ancestor in ancestors[node]:
            children.setdefault(ancestor, set()).add(node)
    return children

"""
Calculate the minimum depth (distance from the closest root) of each term
"""
def get_depth(ontology):
    # Identify nodes with no predecessors
    nodes, parents = get_labels(ontology).keys(), get_parents(ontology)
    roots = set(nodes) - set(parents.keys())
    # Init the dictionary
    depth = {}
    for node in nodes:
        c = 0
        # Get parents of the node, return None if node is a root
        node_parents = parents.get(node)
        while node_parents:
            c += 1
            # Break the loop if the root is among parents
            if roots.intersection(set(node_parents)):
                break
            # Get the parents of current parents (1 level up)
            node_parents = [term for parent in node_parents for term in parents.get(parent, [])]
        depth[node] = c
    return depth

"""
Pipeline for the enrichment test. Take as input two dataframe and the onotlogy file.
Return a Dataframe with as index the GO_ids and values:
1. the p-value and Odd-Ration of the Fisher exact test,
2. the depth computed from the ontology file
3. the description of the GO_ids
"""
def enrich(df1, df2, ontology, col_name_go = 'go_id', col_name_descr='go_descr'):
    # 1. Get Fisher results
    df = fisher_test(df1, df2, col_name_go=col_name_go)
    # 2. Get Depth
    depth = get_depth(ontology)
    # 3. Get description
    labels = get_labels(ontology)
    # 4. Update dataframe
    labels_, depth_ = [], []
    for go_id in df.index:
        labels_.append(labels[go_id])
        depth_.append(depth[go_id])
    df['depth'] = depth_
    df[col_name_descr] = labels_
    # 5. Assign to every GO term the minimum pvalue between its pvalue and its children ones
    df = transmit_pvalue(df, ontology)
    # 6. Return dataframe
    return df

"""
Function that assign to every GO terms the minimum p-value between its own p-value and the p-values of their children.
"""
def transmit_pvalue(enrichment, ontology):
    # 1. Get the children of every GO term
    children_dict = get_children(ontology)
    # 2. For every GO in our enrichment dataset we assign to it the minimum p-value of its children
    for go_id in enrichment.index:
        # Check if the GO term has child
        if children_dict.get(go_id):
            # Retrieve the set of the p-values of all its children
            pvalues = enrichment['p-value'][enrichment.index.isin(children_dict[go_id])]
            # Check we have some children in the dataset. Otherwise we have an empy set 'pvalues'
            if list(pvalues.values):
                # Check if the mimimum pvalue is actually lower than the ancestor one
                min_pvalue = pvalues.min()
                if min_pvalue < enrichment['p-value'][enrichment.index == go_id].values[0]:
                    # If all the conditions are True we assign the minimum pvalue
                    enrichment['p-value'][enrichment.index == go_id] = min_pvalue
    return enrichment

"""
Filter the enrich dataframe by taking out GO_terms with high p-value or high depth
"""
def enrich_filter(df, max_pvalue=0.05, max_depth=5):
    df_filter = df[(df['p-value'] < max_pvalue) & (df['depth'] < max_depth)]
    df_filter['score'] = np.log(1/df['p-value'])
    return df_filter

"""
Create the word cloud of the description of the enriched dataframe, using as frequencies the inverse of p-value
"""
def word_cloud(df, col_name, col_score, *args, **kwargs):
    return WordCloud(*args, **kwargs).generate_from_frequencies({
        row[col_name]: row[col_score] for i, row in df.iterrows()
    })

if __name__ == '__main__':

        # 1. Define arguments
        parser = argparse.ArgumentParser()
        ### Path of the file containing the ontology graph structure
        parser.add_argument('--ontology_path',        type=str,   default='data/go/go.json.gz')
        ### Paths of the two dataframe that must be compared
        parser.add_argument('--original_df_path',     type=str,   default='data/go/go_positive.csv')
        parser.add_argument('--background_df_path',   type=str,   default='data/go/go.csv')
        ### Out directory of the results and of the WordCloud
        parser.add_argument('--out_path',             type=str)
        parser.add_argument('--out_wordcloud',        type=str)
        ### Parameters of the filter
        parser.add_argument('--p_value',              type=float, default=0.05)
        parser.add_argument('--depth',                type=int,   default=4)
        ### Names of the columns of the input dataframe containing the GO_id and the description
        parser.add_argument('--col_name_go_id',       type=str,   default='go_id')
        parser.add_argument('--col_name_descr',       type=str,   default='go_descr')

        # 2. Define dictionary of args
        args = parser.parse_args()

        # 3. Load the required files
        with gzip.open(args.ontology_path) as f:
          ontology = json.load(f)

        ### DF1 and DF2
        original_go= pd.read_table(args.original_df_path,
                                   dtype={'entry_ac': str,
                                          args.col_name_go_id: str,
                                          args.col_name_descr: str})

        background_go= pd.read_table(args.background_df_path,
                                     dtype={'entry_ac': str,
                                            args.col_name_go_id: str,
                                            args.col_name_descr: str})

        # 3. Compute the enrichness
        enrich_result = enrich(df1=original_go,
                               df2=background_go,
                               ontology=ontology,
                               col_name_descr=args.col_name_descr,
                               col_name_go=args.col_name_go_id)

        # 4. Filter the results and create the WordCloud
        ### Results
        enrich_result = enrich_filter(df = enrich_result, max_depth=args.depth, max_pvalue=args.p_value)
        ### WordCloud.
        wc = word_cloud(df=enrich_result, col_name=args.col_name_descr, col_score='score')
        fig = plt.imshow(wc, interpolation='bilinear')


        # 5. Output the results
        ### We save the results if we have a out_dir otherwise we display the results
        if args.out_path:
            enrich_result.to_csv(args.out_path, sep='\t')
        else:
            print(enrich_result)

        ###
        if args.out_wordcloud:
            plt.savefig(fig, args.out_wordcloud)
        else:
            plt.show()
