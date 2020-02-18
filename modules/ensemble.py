####################################
### ENSEMBLE AND MAJORITY VOTING ###
####################################


# Dependencies
import pandas as pd


# Define a method for majority voting
def majority_voting(models_out, threshold=None):
    """
    Input:
        1. models_out:  models output, iterable containing dataframe (entry_ac, positive matching)
        2. threshold:   threshold for majority voting, default number of models
    Output:
        1. set of values passing threshold (predicted positive)
        2. set of values not passing threshold (predicted negative)
    """

    # Set default threshold
    threshold = len(models_out) if threshold is None else threshold

    # Majority voting on accessions
    ensemble_out = pd.concat(models_out)  # Stack all dataframes together
    # Threshold models
    entry_ac = ensemble_out['entry_ac']
    counts = entry_ac.value_counts()  # Get counts for each unique value
    # Apply threshold, restrieve binary mask and apply it to entries list
    passed = counts.index[counts >= threshold].tolist()
    # Save only models which have passed the given threshold
    ensemble_out = ensemble_out[entry_ac.isin(passed)]

    # Create a list of set for every grouped accession
    ensemble_out = ensemble_out.groupby(by='entry_ac').agg({
        'positive': lambda x: list(x)
    }).reset_index()
    # Turn list of sets into list
    ensemble_out['positive'] = ensemble_out['positive'].apply(lambda matches: [pos for match in matches for pos in match])
    # Apply threshold on list values
    ensemble_out['positive'] = ensemble_out['positive'].apply(lambda pos: set(p for p in set(pos) if pos.count(p) >= threshold))

    # Return majority voting output
    return ensemble_out
