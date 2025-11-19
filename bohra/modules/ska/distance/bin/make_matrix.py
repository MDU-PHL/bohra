#!/usr/bin/env python3

import pandas as pd
import sys

def make_matrix(df, col1, col2, vals) -> pd.DataFrame:
    """
    Create a matrix from two columns of a DataFrame and a list of values.

    Args:
        df (pd.DataFrame): The input DataFrame.
        col1 (str): The name of the first column.
        col2 (str): The name of the second column.
        vals (list): A list of values to include in the matrix.

    Returns:
        pd.DataFrame: A DataFrame representing the matrix.
    """
    idx = sorted(set(df[col1]).union(df[col2]))
    # Create a pivot table
    matrix = df.pivot(index = col1, columns = col2, values = vals).reindex(index = idx, columns = idx).fillna(0, downcast = 'infer').pipe(lambda x: x+x.values.T)
    # Reindex to include all values
    matrix = matrix.reset_index()
    matrix = matrix.rename(columns = {col1: 'Isolate'})

    return matrix


pth = sys.argv[1]
df = pd.read_csv(pth, sep = '\t')
matrix = make_matrix(df, 'Sample1', 'Sample2', 'Distance')
matrix.to_csv("distances.tsv", sep = '\t', index = False)