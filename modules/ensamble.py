import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import parser_lib
from itertools import product

"""
Get the union of the passed proteins of the HMM and PSSM (Get proteins that passed PSSM or HMM).
"""
def ensamble(pssm_path, hmm_path, true_positive_path):
    # Get PSSM results
    df_pssm = pd.read_csv(pssm_path, index_col = [0])
    index_pssm = set(df_pssm.index)
    # Get HMM results
    df_hmm = pd.read_csv(hmm_path, index_col = [0, 1])
    index_hmm = set([index[0] for index in df_hmm.index])
    # Get true_positive
    true_positive_prot = parser_lib.fasta_to_dict(path=true_positive_path)
    true_positive = set([key.split('|')[1]for key in true_positive_prot.keys()])


    # Get the union of the two results
    ensamble_positive = index_pssm.union(index_hmm)
    ensamble_df = pd.DataFrame()
    PSSM = [x in index_pssm for x in ensamble_positive]
    HMM = [x in index_hmm for x in ensamble_positive]
    TP = [x in true_positive for x in ensamble_positive]
    ensamble_df['PSSM'] = PSSM
    ensamble_df['HMM'] = HMM
    ensamble_df['TP'] = TP
    ensamble_df.index = ensamble_positive
    return ensamble_df

"""
Create the confusion matrix and compute some statistics (precision, recall, accuracy)
using the results of the ensamble and the positive and negative test set
(e.g results/positive_pfam.fasta results/negative_pfam.fasta).
Return a numpy matrix for the confusion_matrix and dictionary
"""
def confusion_matrix(ensamble_df, positve_path, negative_path):
    # Positive IDs Sets
    positive_prot = parser_lib.fasta_to_dict(path=positve_path)
    positive = set([key.split('|')[1]for key in positive_prot.keys()])
    # Negative IDs Sets
    negative_prot = parser_lib.fasta_to_dict(path=negative_path)
    negative = set([key.split('|')[1]for key in negative_prot.keys()])
    # True Positive and False Positive IDs Sets
    true_positive = set(df.index[df.TP == True])
    false_positive = set(df.index[df.TP == False])
    # True Negative and False Negative IDs Sets
    true_negative = negative.difference(false_positive)
    false_negative = positive.difference(true_positive)
    # True and False Positive and True and False Negative
    tp, fp, tn, fn = len(true_positive), len(false_positive), len(true_negative), len(false_negative)
    # Confusion Matrix
    cf = np.matrix([[tp, fp], [fn, tn]])
    # Statistics
    stats = {}
    stats['precision'] = tp/(tp+fp)
    stats['recall'] = tp/(tp+fn)
    stats['accuracy'] = (tp+tn)/(tp+fp+fn+tn)
    stats['weighted_accuracy'] = (tp/(tp+fn) + tn/(tn+fp))/2
    return cf, stats


if __name__ == '__main__':
    # 1. Define the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--pssm_df_path',     type=str,   default='results/PSSM_df.csv')
    parser.add_argument('--hmm_df_path',      type=str,   default='results/HMM_df.csv')
    parser.add_argument('--positive_path',    type=str,   default='data/positive_pfam.fasta')
    parser.add_argument('--negative_path',    type=str,   default='data/negative_pfam.fasta')
    parser.add_argument('--out_path',         type=str,   default='results/ensamble/ensamble_df.csv')
    parser.add_argument('--out_path_img',     type=str,   default='results/ensamble/confusion_matrix.png')
    parser.add_argument('--out_path_metrics', type=str,   default='results/ensamble/metrics.csv')
    parser.add_argument('--verbose',          type=bool,  default= True)

    # 2. Generate dictionary of args
    args = parser.parse_args()

    # 3. Run the ensemble
    df = ensamble(args.pssm_df_path, args.hmm_df_path, args.positive_path)
    df.to_csv(args.out_path)

    # 4. Generate statistics and confusion matrix
    cf, stats = confusion_matrix(df, args.positive_path, args.negative_path)

    # 5. Save metrics
    metrics = pd.Series([stats[key] for key in stats], index=stats.keys())
    metrics.to_csv(args.out_path_metrics)

    # 6. Show confusion matrix
    cf_log = np.log(cf) # We acutally use the log because the differences are too high
    plt.imshow(cf_log, interpolation='nearest', cmap=plt.cm.Blues)
    # Set colorbar
    plt.colorbar(plt.cm.ScalarMappable(cmap=plt.cm.Blues),
                 boundaries = [cf[1, 1]*i/30 for i in range(33)],
                 fraction=0.04, pad=0.1)
    # Set ticks
    tick_marks = range(2)
    plt.xticks(tick_marks, ['True', 'False'], rotation=0)
    plt.yticks(tick_marks, ['True', 'False'])
    # Print text
    for i, j in product(range(cf.shape[0]), range(cf.shape[1])):
        plt.text(j, i, format(cf[i, j], 'd'),
                 horizontalalignment="center", fontsize=18)
    # Set labels
    plt.title('Confusion Matrix')
    plt.ylabel("Predicted\n--\nModel Hit")
    plt.xlabel("PFAM: 00397\n--\nReal")
    # Save and show
    plt.savefig(args.out_path_img)
    plt.show()

    # 7. Print results if option verbose is True
    if args.verbose:
        print('Dataframe of hits:\n{}\n'.format(df))
        print('Confusion Matrix:\n{}\n'.format(cf))
        print('Metrics:\n{}'.format(metrics))
