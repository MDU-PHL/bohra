#!/usr/bin/env python3

import pandas as pd
import sys



def combine_stats(dist_df, stats_df) -> pd.DataFrame:
    dist_df = dist_df.rename(columns = {"Match count":"# shared kmers"})
    # print(dist_df)
    tmp1 = stats_df.rename(columns = {'Isolate': 'Sample1', '# kmers': 'Sample1 #kmers'})
    # print(tmp1)
    tmp2 = stats_df.rename(columns = {'Isolate': 'Sample2', '# kmers': 'Sample2 #kmers'})
    # print(tmp2)
    merged = dist_df.merge(tmp1[['Sample1', 'Sample1 #kmers']], on = 'Sample1', how = 'left')
    # print(merged)
    merged = merged.merge(tmp2[['Sample2', 'Sample2 #kmers']], on = 'Sample2', how = 'left')
    # print(merged)
    merged["Sample1 % contribution"] = round(merged['# shared kmers'] / merged['Sample1 #kmers'] * 100, 2)
    merged["Sample2 % contribution"] = round(merged['# shared kmers'] / merged['Sample2 #kmers'] * 100, 2)
    merged["Distance"] = merged["Distance"].apply(lambda x: int(x))
    # print(merged)
    merged["Minimum % contribution"] = merged[["Sample1 % contribution", "Sample2 % contribution"]].min(axis=1)
    merged = merged[[ "Sample1", "Sample1 #kmers", "Sample2", "Sample2 #kmers", "# shared kmers", "Sample1 % contribution", "Sample2 % contribution", "Minimum % contribution","Distance"]]
    print(merged)
    merged.to_csv("ska_distance_stats.tsv", sep = '\t', index = False)
    return merged
def open_stats(ska_nk):

    try:

        stats = pd.read_csv(ska_nk, sep = '\t')
        return stats
    except Exception as e:
        print(f"Error reading SKA nk stats file: {e}")
        sys.exit(1)
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


dist_stats = sys.argv[1]
df = pd.read_csv(dist_stats, sep = '\t')
ska_nk = sys.argv[2]

stats = open_stats(ska_nk)
distance_metrics = combine_stats(df, stats)
matrix = make_matrix(df, 'Sample1', 'Sample2', 'Distance')
matrix.to_csv("distances.tsv", sep = '\t', index = False)