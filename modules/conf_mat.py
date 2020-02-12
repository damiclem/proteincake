##################################
### RESULTS EVALUATION LIBRARY ###
##################################


# Dependencies
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


# Generate confusion matrix
def compute(positive, negative, pred_positive, pred_negative):
# def get_conf_matrix(ensemble_df, positve_path, negative_path):
    """
    Create the confusion matrix and compute some statistics (precision, recall, accuracy)
    Return a numpy matrix for the confusion_matrix along with computed statistics
    Input:
        1. positive:        set of actually positive values
        2. negative:        set of actually negative values
        3. pred_positive:   set of values predicted positive
        4. pred_negative:   set of values predicted negative
    Output:
        1. conf_mat: confusion matrix (numpy ndarray with shape (2, 2))
        2. precision:
        3. recall:
        4. accuracy:
    """
    # Compute confusion matrix entries
    tp = len(positive & pred_positive)
    fp = len(negative & pred_positive)
    tn = len(negative & pred_negative)
    fn = len(positive & pred_negative)
    # Fill confusion matrix entries
    conf_mat = np.empty(shape=(2,2), dtype=np.int)
    conf_mat[1, 0], conf_mat[1, 1] = tp, fp
    conf_mat[0, 0], conf_mat[0, 1] = fn, tn

    # Compute statistics
    precision = tp/(tp+fp)
    recall = tp/(tp+fn)
    accuracy = (tp+tn)/(tp+fp+fn+tn)
    weighted_accuracy = (tp/(tp+fn) + tn/(tn+fp))/2

    # Return either the matrix and the statistics
    return conf_mat, precision, recall, accuracy, weighted_accuracy


# Show confusion matrix
def plot(positive, negative, pred_positive, pred_negative, ax=None):

    # Get confusion matrix, along with statistics
    conf_mat, precision, recall, accuracy, weighted_accuracy = compute(
        positive, negative, pred_positive, pred_negative
    )

    # Format output text
    out_mat = conf_mat.astype(np.unicode_)
    for i in range(out_mat.shape[0]):
        for j in range(out_mat.shape[1]):
            out_mat[i, j] = format(conf_mat[i, j], 'd')

    # Create heatmap
    _ = sns.heatmap(np.log10(1 + conf_mat), cbar=False, annot=out_mat,
                    cmap=plt.cm.Blues, fmt='',
                    annot_kws={'size': 16, 'ha': 'center', 'va': 'center'},
                    ax=ax)
    # Set axis
    _ = ax.set_ylim(0, conf_mat.shape[0])
    # Set ticks
    _ = ax.set_xticklabels(['True', 'False'])
    _ = ax.set_yticklabels(['False', 'True'])
    # Add titles
    _ = ax.set_xlabel('Real values')
    _ = ax.set_ylabel('Predicted values')
