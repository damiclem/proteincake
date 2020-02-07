import pandas as pd

"""
Parse the HMM results such that we return a Dataframe with multindex (ID, domain number) and values some statistics about
the domain.
"""
def hmm_parser(path):
    with open(path) as results:
        # Analyze HMM model
        # Columns
        results.__next__()
        col = results.__next__().split()[14:-3]
        col[4], col[5] = 'hmm_{}'.format(col[4]), 'hmm_{}'.format(col[5])
        col[6], col[7] = 'ali_{}'.format(col[6]), 'ali_{}'.format(col[7])
        col[8], col[9] = 'env_{}'.format(col[8]), 'env_{}'.format(col[9])
        # Init MultiIndex and Values
        outer, inner = [], []
        values = []
        # We cycle in the rows to get ID, domain and values
        for lines in results:
            if lines[0:2] == 'sp':
                line = lines.split()
                # Here we get the values
                values.append(line[11:22])
                # Here we get the multindex
                key = line[0].split('|')[1]
                if key not in outer:
                    n = 1
                    outer.append(key)
                    inner.append(1)
                else:
                    outer.append(key)
                    inner.append(n)
                n += 1
        # Create the MultiIndex
        arrays = [outer, inner]
        index = pd.MultiIndex.from_arrays(arrays, names=('ID', 'Domain'))
        # Return df
        df = pd.DataFrame(values, columns=col, index=index)
        # Transform values from str to float
        df.iloc[:, :4] = df.iloc[:, :4].applymap(lambda x: float(x))
        return df



"""
Parse the PSSM results such that we return a Dataframe with index ID and Bit-Score and E-Value.
"""
def pssm_parser(path):
    with open(path) as results:
        # Analyse pssm model
        while results.__next__()[0:9] != 'Sequences':
            continue
        df = pd.DataFrame()

        condition = True
        results.__next__()
        while(condition):
            line = results.__next__().split()
            if len(line) > 0:
                df[line[0]] = [line[-2], line[-1]]
            else:
                condition = False
        df = df.transpose()
        df.columns = ['Bit-Score', 'E-Value']
        return df

"""
Take a path of a fasta file and return a dict {id: sequence} (tuned on the example of uniref90)
"""
def fasta_to_dict(path):
    with open(path, 'r') as file:
        lines = [line for line in file]

    # Get two separate list
    ids = []
    chain = []
    temp=[]
    for n, l in enumerate(lines):
        if l[0] == '>':
            sentence = l.split()[0]
            ids.append(sentence[1:]) # We take only the first part as identifier
            if n>1:
                chain.append(temp)
                temp = []
        else:
            temp.append(l[:-1]) # We remove the '\n' at the end
    chain.append(temp)
    # Get a dictionary
    dict_fasta = {key: value for key, value in zip(ids, chain)}
    for key in dict_fasta:
        temp = ''
        n = len(dict_fasta[key])
        for i in range(n):
            temp += dict_fasta[key][i]
        dict_fasta[key] = temp

    return dict_fasta