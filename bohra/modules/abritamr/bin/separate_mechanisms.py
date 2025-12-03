#!/usr/bin/env python3
import sys,  json
import pandas as pd


def get_snps(col:str) -> str:
    snps = []
    for mech in col.split(','):
        if '_' in mech:
            snps.append(mech)
            # return mech
    return ','.join(snps)


def get_genes(col:str) -> str:
    genes = []
    for mech in col.split(','):
        if '_' not in mech:
            genes.append(mech)
            # return mech
    return ','.join(genes)



df = pd.read_csv(sys.argv[1], sep='\t')
replace_negs = {
    "CPase_ESBL_AmpC_16S_NEG":"",
    "Mec_VanAB_Linez_NEG":"",
    "Van_Linez_NEG":"",
    "CPase_16S_mcr_NEG":""
}
df = df.replace(replace_negs)
df['Reportable AMR genes'] = df['Reportable AMR mechanisms'].apply(lambda x: get_genes(x))
df['Reportable AMR SNPs'] = df['Reportable AMR mechanisms'].apply(lambda x: get_snps(x))
df['Other AMR genes'] = df['Other AMR mechanisms'].apply(lambda x: get_genes(x))
df['Other AMR SNPs'] = df['Other AMR mechanisms'].apply(lambda x: get_snps(x))

df[['Isolate',
    'Reportable AMR genes',
    'Reportable AMR SNPs',
    'Other AMR genes',
    'Other AMR SNPs',
    'Species']].to_csv("reportable_amr_matches.txt", sep='\t', index=False)