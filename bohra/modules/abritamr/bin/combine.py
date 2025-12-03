#!/usr/bin/env python3
import sys,  json
import pandas as pd


def get_genes(_file, _type, _dict):
    
    with open(_file, 'r') as f:
        data = f.read().strip().split('\n')
        header = data[0].split('\t')
        line = data[1].split('\t')
        for i in range(1,len(header)):
            # print(header[i])
            l = [gn for gn in line[i].split(',')] if _type != 'partials' else [f"{gn}^" for gn in line[i].split(',')]
            # l = line[i] if _type != 'partials' else f"{line[i]}^"
            # print(line[i])
            # print(l)
            if header[i] in _dict:
                    _dict[header[i]].extend(l)
            else:
                _dict[header[i]] = l
        print(_dict)
        return _dict

def extract_contigs_plasmids(plasmid: pd.DataFrame) -> list:
    
    contigs =[ f"{i.split()[0]}" for i in plasmid[plasmid["molecule_type"] == "plasmid"]['contig_id'].unique()]
    return contigs

def extract_genes_on_plasmid(amrfinder: pd.DataFrame, contigs: list) -> list:
    
    genes_on_plasmid = amrfinder[amrfinder["Contig id"].isin(contigs)]["Gene symbol"].unique()
    return genes_on_plasmid

def get_plasmid_designation_species(plasmid: pd.DataFrame, ctg:str) -> tuple:

    plasmid_designation = plasmid[plasmid["contig_id"].str.contains(ctg)]["mash_nearest_neighbor"].unique()[0]
    plasmid_species = plasmid[plasmid["contig_id"].str.contains(ctg)]["mash_neighbor_identification"].unique()[0]
    return plasmid_designation, plasmid_species

def make_amr_df(genes_dict: dict, genes_on_plasmid:list, plasmid:pd.DataFrame, amrout:pd.DataFrame) -> pd.DataFrame:
    
    for row in genes_dict:
        print(row)
        
        if row != 'Isolate':
            print(genes_dict[row])

            res = []
            for gene in genes_dict[row]:
                if gene in genes_on_plasmid:
                    ctg = amrout[amrout["Gene symbol"] == gene]["Contig id"].unique()[0]
                    plasmid_designation, plasmid_species = get_plasmid_designation_species(plasmid, ctg)
                    gene = f"{gene}[P]"
                
                res.append(gene)
            genes_dict[row] = ",".join(res)
    print(pd.DataFrame(genes_dict, index = [0]))
    return pd.DataFrame(genes_dict, index = [0])

def combine_plasmid_data(plasmid: pd.DataFrame, amrfinder: pd.DataFrame, contigs: list, sid:str) -> pd.DataFrame:
    
    data = {"Isolate": sid,
            "_children" : [],
            "mash_nearest_neighbor":"",
            "mash_neighbor_identification":"",
            "Present_on_contigs":"",
            "AMR_genes":"",
            "rep_type":""}
    plasmid = plasmid[plasmid["molecule_type"] == "plasmid"]
    
    mnm = []
    mni = []
    poc = []
    ag = []
    rt = []


    for p in plasmid["mash_nearest_neighbor"].unique():
        mnm.append(p)
       
        tmp = plasmid[plasmid["mash_nearest_neighbor"] == p]
        ctgs = [f"{i.split()[0]}" for i in tmp["contig_id"].unique()]
        mni.append(tmp["mash_neighbor_identification"].unique()[0])
        poc.append(",".join(ctgs))
        ag.append(",".join(sorted(amrfinder[amrfinder["Contig id"].isin(ctgs)]["Gene symbol"].unique())))
        rt.append(",".join(tmp["rep_type(s)"].unique()))
        child = {"Isolate": sid,}
        child["Plasmid accession"] = p
        child["Plasmid species"] = tmp["mash_neighbor_identification"].unique()[0]
        child["Present on contigs"] = ",".join([f"{i.split()[0]}" for i in tmp["contig_id"].unique()])
        child["Replicon type"] = ",".join(tmp["rep_type(s)"].unique())
        child["AMR genes"] = ",".join(sorted(amrfinder[amrfinder["Contig id"].isin(ctgs)]["Gene symbol"].unique()))
        data["_children"].append(child)
    data["Plasmid accession"] = "|".join(mnm) if mnm != [] else ""
    data["Plasmid species"] = "|".join(mni) if mni != [] else ""
    data["Present on contigs"] = "|".join(poc) if poc != [] else ""
    data["AMR genes"] = "|".join(ag) if ag != [] else ""
    data["Replicon type"] = "|".join(rt) if rt != [] else ""
    
    with open(f"plasmid.json", 'w') as f:
        json.dump(data, f, indent = 4)
    # return plasmid

_dict = {'Isolate': sys.argv[1]}

# matches
_dict = get_genes(_file = sys.argv[2], _type = 'matches', _dict = _dict)
# partials 
_dict = get_genes(_file = sys.argv[3], _type = 'partials', _dict = _dict)

amrfinder = pd.read_csv(sys.argv[4], sep = '\t')
plasmid = pd.read_csv(sys.argv[5], sep = '\t')

genes_on_plasmid = extract_genes_on_plasmid(amrfinder, extract_contigs_plasmids(plasmid))
contigs = extract_contigs_plasmids(plasmid)

combine_plasmid_data(plasmid, amrfinder, contigs, _dict['Isolate'])

df_amr = make_amr_df(_dict, genes_on_plasmid, plasmid, amrfinder)
df_amr.to_csv("resistome.txt", sep = '\t', index = False)