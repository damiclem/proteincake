####################################
### ENSEMBLE AND MAJORITY VOTING ###
####################################


# Dependencies
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.cluster import DBSCAN


# Define a method for majority voting
def majority_voting(models_out, threshold=None):
    """
    Input:
        1. models_out:      models output, iterable containing dataframe (entry_ac, positive matching)
        2. threshold:       threshold for majority voting, default number of models
    Output:
        1. ensemble_out:    set of values passing threshold (predicted positive)
    """
    # # Set default threshold
    # threshold = len(models_out) if threshold is None else threshold

    # Concatenate models output
    models_out = pd.concat(models_out, ignore_index=True)  # Stack all models output together
    # Loop thorugh each accession
    for accession in set(models_out.entry_ac.unique().tolist()):
        # Get models with respect to current accession
        models_curr = models_out[models_out.entry_ac == accession]
        # Get indices of currently sliced models
        index_curr = models_curr.index
        # Compute distance matrix as accuracy matrix (symmetric)
        # NOTE: only upper triangular matrix is computed
        nrows = ncols = models_curr.shape[0]  # Define number of rows and columns
        # Define clusters container
        labels = 0
        # Handle case single cell matrix
        if nrows > 1:
            # Initialize distance matrix
            distance_matrix = np.zeros(shape=(nrows, ncols), dtype=np.float)
            # Fill distance matrix
            for i in range(nrows):
                for j in range(i + 1, ncols):
                    # Define "truncated accuracy" for i-th row and j-th column
                    tp = len(models_curr.loc[index_curr[i], 'positive'] & models_curr.loc[index_curr[j], 'positive'])  # Intersection
                    f = len(models_curr.loc[index_curr[i], 'positive'] | models_curr.loc[index_curr[j], 'positive'])  # Union
                    # Fill distance matrix
                    distance_matrix[i, j] = distance_matrix[j, i] = 1 - (tp / f)
            # Compute cluster labels through DBScan
            labels = DBSCAN(eps=0.4, min_samples=1).fit_predict(distance_matrix)
        # Assign model output to clusters
        models_out.loc[index_curr, 'label'] = labels

    # For each label in each accession, compute majority voting
    # Group each cluster and merge together positive values
    ensemble_out = models_out.groupby(by=['entry_ac', 'label']).agg({
        'positive': {
            'positive': lambda matches: [pos for match in list(matches) for pos in match],
            'count': 'count'
        }
    }).reset_index()
    # Reset column names
    ensemble_out.columns = ['entry_ac', 'label', 'positive', 'count']
    # Compute majority voting
    ensemble_out['positive'] = ensemble_out.apply(
        lambda x: set(p for p in set(x['positive']) if x['positive'].count(p) >= (x['count'] / 2)),
        axis=1
    )

    # Get only rows with at least matching position in domain
    ensemble_out = ensemble_out[ensemble_out['positive'].apply(len) > 0]
    # Slice columns
    ensemble_out = ensemble_out[['entry_ac', 'positive']]
    # return result
    return ensemble_out


# Get ensemble output from multiple predictions
if __name__ == '__main__':

    # 1. Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--models_out', type=str, nargs='+', required=True)
    parser.add_argument('--threshold', type=int, required=False)
    parser.add_argument('--out_path', type=str)
    args = parser.parse_args()

    # 2. Load models from csv
    models_out = [pd.read_csv(model_out, sep='\t') for model_out in models_out]

    # 3. Run majority voting
    ensemble_out = majority_voting(models_out, args.threshold)

    # 4. Handle result
    if args.out_path:
        # Case output path has been defined: store to file
        ensemble_out.to_csv(args.out_path, index=False)
    else:
        # Case no output path defined: print out
        print(ensemble_out.to_string())
