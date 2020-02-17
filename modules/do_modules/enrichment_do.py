import json
import gzip
import argparse
import warnings
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from tqdm import tqdm
from scipy.stats import fisher_exact
from wordcloud import WordCloud

#warnings.filterwarnings("ignore")
pd.options.mode.chained_assignment = None



"""
To retrieve the GO that are parents, we cycle over ontology["graphs"][0]["edges"] which is a list of dictionary.
Every dictionary is about a GO id (sub) with a relation (pred) with another GO (obj).
We create a dictionary (parents) with as keys the sons and as values the parents.
"""
def get_parents(ontology):
    return ontology.is_a[ontology.is_a.isna() == False].to_dict()

"""
Here we cycle over the nodes to obtain a dictionary of GO_id with as value a description.
- ontology["graphs"][0]["nodes"] is a list of dictionary with dict_keys(['id', 'meta', 'type', 'lbl'])
- ontology["graphs"][0]["nodes"][1]['lbl'] is the value (e.g: "endocytosed synaptic vesicle processing via endosome")
"""
def get_labels(ontology):
    return ontology['name'].to_dict()

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
    nodes, parents = ontology.do_id, get_parents(ontology)
    roots = set(nodes) - set(parents.keys())
    # Init the dictionary
    depth = {}
    for node in tqdm(nodes, ncols=100,
                     bar_format='{l_bar}{bar:40}{r_bar}{bar:-40b}',
                     desc='Depth            '):
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
Perform Fisher test. An Odd-Ratio above 77 tells us the GO prefers the first dataframe (p-value < 0.05),
while an Odd-Ratio under 0.013 tells us the GO prefers the second dataframe.
Return a Dataframe with index the GO and values the Odd-Ratio and the p-value.
"""
def fisher_test(df1, df2, col_name_do = 'do_id'):

    # Inint dict
    results = {}

    # Get the number of occurrances of the GO counts
    dict1, dict2 = dict(df1[col_name_do].value_counts()), dict(df2[col_name_do].value_counts())

    # Compute the intersaction of the GO terms
    key_intersection = set(dict1.keys()).intersection(set(dict2.keys()))

    for key in tqdm(key_intersection, ncols=100,
                    bar_format='{l_bar}{bar:40}{r_bar}{bar:-40b}',
                    desc='Fisher Test      '):
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
Function that assign to every GO terms the minimum p-value between its own p-value and the p-values of their children.
"""
def transmit_pvalue(enrichment, ontology):
    # 1. Get the children of every GO term
    children_dict = get_children(ontology)
    # 2. For every GO in our enrichment dataset we assign to it the minimum p-value of its children
    for do_id in tqdm(enrichment.index, ncols=100,
                      bar_format='{l_bar}{bar:40}{r_bar}{bar:-40b}',
                      desc='Propagate p-value'):
        # Check if the GO term has child
        if children_dict.get(do_id):
            # Retrieve the set of the p-values of all its children
            pvalues = enrichment['p-value'][enrichment.index.isin(children_dict[do_id])]
            # Check we have some children in the dataset. Otherwise we have an empy set 'pvalues'
            if list(pvalues.values):
                # Check if the mimimum pvalue is actually lower than the ancestor one
                min_pvalue = pvalues.min()
                if min_pvalue < enrichment['p-value'][enrichment.index == do_id].values[0]:
                    # If all the conditions are True we assign the minimum pvalue
                    enrichment['p-value'][enrichment.index == do_id] = min_pvalue
    return enrichment



"""
Pipeline for the enrichment test. Take as input two dataframe and the onotlogy file.
Return a Dataframe with as index the GO_ids and values:
1. the p-value and Odd-Ration of the Fisher exact test,
2. the depth computed from the ontology file
3. the description of the GO_ids
"""
def enrich(df1, df2, ontology, col_name_do = 'do_id', col_name_descr='name'):
    # 1. Get Fisher results
    df = fisher_test(df1, df2, col_name_do=col_name_do)
    # 2. Get Depth
    depth = get_depth(ontology)
    # 4. Update dataframe
    labels_, depth_ , do_found= [], [], []
    for do_id in df.index:
         if depth.get(do_id):
            do_found.append(do_id)
            labels_.append(ontology[col_name_descr][ontology[col_name_do] == do_id].values[0])
            depth_.append(depth[do_id])

    df = df[df.index.isin(do_found)]
    df['depth'] = depth_
    df[col_name_descr] = labels_
    df = transmit_pvalue(df, ontology)
    # 5. Return dataframe
    return df

"""
Filter the enrich dataframe by taking out GO_terms with high p-value or high depth
"""
def enrich_filter(df, max_pvalue=0.05, max_depth=5):
    df_filter = df[(df['p-value'] < max_pvalue) & (df['depth'] <= max_depth)]
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
        parser.add_argument('--ontology_path',        type=str,   default='data/do/do_ontology.csv')
        ### Paths of the two dataframe that must be compared
        parser.add_argument('--target_path',          type=str,   default='data/do/do_original.csv')
        parser.add_argument('--background_path',      type=str,   default='data/do/do.csv')
        ### Out directory of the results and of the WordCloud
        parser.add_argument('--out_path',             type=str)
        parser.add_argument('--out_wordcloud',        type=str)
        ### Parameters of the filter
        parser.add_argument('--p_value',              type=float, default=0.05)
        parser.add_argument('--depth',                type=int,   default=4)
        parser.add_argument('--bonferroni',           type=int,  default=1)
        ### Names of the columns of the input dataframe containing the GO_id and the description
        parser.add_argument('--col_name_do_id',       type=str,   default='do_id')
        parser.add_argument('--col_name_descr',       type=str,   default='name')

        # 2. Define dictionary of args
        args = parser.parse_args()
        # 3. Load the required files
        ontology = pd.read_csv(args.ontology_path, dtype=str, sep='\t', index_col=[0])
        ontology.index = ontology.do_id.values
        ### DF1 and DF2
        original_do= pd.read_table(args.target_path, dtype=str, index_col = [0])
        background_do= pd.read_table(args.background_path, dtype=str, index_col = [0])

        # 3. Compute the enrichness
        enrich_result = enrich(df1=original_do,
                               df2=background_do,
                               ontology=ontology,
                               col_name_descr=args.col_name_descr,
                               col_name_do=args.col_name_do_id)

        # Use bonferroni
        if args.bonferroni:
            if ((enrich_result['p-value'] > args.p_value/enrich_result.shape[0]) | (enrich_result['depth'] > args.depth)).all():
                warnings.warn('No object passed the filter. Returning non-filtered dataset.\nClosing...')
                enrich_result['score'] = np.zeros((1, enrich_result.shape[0])).reshape(-1) - 1
                if args.out_path:
                    enrich_result.to_csv(args.out_path, sep='\t')
                else:
                    print(enrich_result)
                sys.exit(1)


            # 4. Filter the results and create the WordCloud
            ### Results
            enrich_result = enrich_filter(df = enrich_result, max_depth=args.depth, max_pvalue=args.p_value/enrich_result.shape[0])
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
                plt.savefig(args.out_wordcloud)
            else:
                plt.show()

        # Don't use bonferroni
        else:
            if ((enrich_result['p-value'] > args.p_value) | (enrich_result['depth'] > args.depth)).all():
                warnings.warn('No object passed the filter. Returning non-filtered dataset.\nClosing...')
                enrich_result['score'] = np.zeros((1, enrich_result.shape[0])).reshape(-1) - 1
                if args.out_path:
                    enrich_result.to_csv(args.out_path, sep='\t')
                else:
                    print(enrich_result)
                sys.exit(1)


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
                plt.savefig(args.out_wordcloud)
            else:
                plt.show()
