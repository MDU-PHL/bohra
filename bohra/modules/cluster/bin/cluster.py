#!/usr/bin/env python3
import numpy,sys, pathlib, pandas
from sklearn.cluster import AgglomerativeClustering



def extract_unclustered(df, col, col1):

    u = df.groupby(col)[col1].agg('count').reset_index()

    uc = list(u[u[col1] == 1][col])

    return uc


thresholds = [int(i) for i in sys.argv[2].split(',')]
# get distances into lists
tabfile = sys.argv[1]
if pathlib.Path(tabfile).exists():
    isos = []
    mat = []
    if isinstance(tabfile, str):
        print(f"opening {tabfile}")
        with open(tabfile, 'r') as d:
            tab = d.read().split('\n')

            # get isolates
            for line in tab[1:len(tab)]:
                ln = line.split('\t')
                if len(ln) > 1:
                    # isos.append(ln[0])
                    mat.append(ln[1:len(ln)])

            isos = tab[0].split('\t')[1:]
    print(len(mat))
    clusters = pandas.DataFrame()
    # convery the lists to a numpy array
    print("converting matrix to numpy array")
    X = numpy.array(mat)
    X = X.astype(numpy.float64)

    for level in thresholds:
        print(f"Clustering at {level}")
        clustering = AgglomerativeClustering(n_clusters = None, metric = 'precomputed',linkage = f"{sys.argv[3]}", distance_threshold =int(level)).fit(X)
        df = pandas.DataFrame(data = {'ID': isos, f"Tx:{level}": clustering.labels_})
        print(len(clustering.labels_))
        print(len(isos))
        print(df.shape) 
        colname = f"Tx:{level-1}"
        df = df.rename(columns={f"Tx:{level}": colname})
        print(len(df[f"{colname}"].unique()))
        df[f"{colname}"] = df[f"{colname}"] + 1
        df = df.fillna('')
        df[f"{colname}"] = df[f"{colname}"].apply(lambda x: f"{x}")
        print(df.shape)

        uc = extract_unclustered(df = df, col = f"{colname}", col1 = 'ID')
        df[f"{colname}"] = numpy.where(df[f"{colname}"].isin(uc), 'UC', df[f"{colname}"])
        df = df.drop_duplicates(subset=['ID'])
        if clusters.empty:
            clusters = df
        else:
            clusters = df.merge(clusters, right_on = 'ID',left_on ='ID', how = 'inner')
        clusters = clusters.drop_duplicates()
        print(clusters.shape)
    clusters = clusters.fillna('')

    clusters.to_csv(f'clusters.txt', sep = '\t', index = False)
