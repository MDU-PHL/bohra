#!/usr/bin/env python3
import numpy as np
import pandas as pd
import sys

def get_lower_triangular_matrix(matrix, isos):
    """
    Converts a square matrix to a lower triangular matrix.

    Args:
        matrix: A list of lists representing the square matrix.

    Returns:
        A list of lists representing the lower triangular matrix.
    """
    
    rows = len(matrix)

    lower_triangle = [[len(isos)], [isos[0]]]
    for i in range(1,rows):
    
        r = [isos[i]]
        for j in range(i + 1):
            
            if i != j:
            # print(matrix[i][j])
                r.append(matrix[i][j])
        lower_triangle.append(r)

    
    lines = []
    for row in lower_triangle:
        line = '\t'.join(map(str, row))
        lines.append(line)
    


    return '\n'.join(lines)
        #     lower_triangular[i][j] = matrix[i][j]
    # return lower_triangular

#Example
matrix = pd.read_csv(sys.argv[1], sep='\t')
isos = matrix.iloc[:, 0].values
iso_col = matrix.columns[0]
lower_triangle = get_lower_triangular_matrix(matrix.set_index(iso_col).values.tolist(), isos)
print(lower_triangle)

