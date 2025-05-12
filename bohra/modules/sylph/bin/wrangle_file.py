#!/usr/bin/env python3

import sys,pandas as pd

sid = sys.argv[1]
sylph = sys.argv[2]

df = pd.read_csv(sylph, sep="\t")

df["Isolate"] = sid
df["Species"] = df["Contig_name"].apply(lambda x: ' '.join(x.split(" ")[1:3]))
df[["Isolate", "Species", "Taxonomic_abundance", "Sequence_abundance", "Adjusted_ANI", "Naive_ANI"]]\
    .to_csv(sys.stdout, sep="\t", index=False)