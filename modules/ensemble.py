####################################
### ENSEMBLE AND MAJORITY VOTING ###
####################################


# Define a method for majority voting
def majority_voting(models_out, threshold=None):
    """
    Input:
        1. models_out:  models output, iterable containing tuples (positive, negative)
        2. threshold:   threshold for majority voting, default number of models
    Output:
        1. set of values passing threshold (predicted positive)
        2. set of values not passing threshold (predicted negative)
    """

    # Set default threshold
    threshold = len(models_out) if threshold is None else threshold

    # Create a list of values
    pp = [positive for out in models_out for positive in out[0]]  # Predicted positive
    pn = [negative for out in models_out for negative in out[1]]  # Predicted negative

    # Count poisitive values
    pp_counts = dict()
    for positive in pp:
        # Fill dictionary entry of predicted positive counts
        pp_counts[positive] = pp_counts.get(positive, 0) + 1

    # Filter positive values using threshold
    pp = set([k for k, c in pp_counts.items() if c >= threshold])
    pn = set([k for k, c in pp_counts.items() if c < threshold]) | set(pn)

    # Return predicted positive and predicted negative
    return pp, pn
