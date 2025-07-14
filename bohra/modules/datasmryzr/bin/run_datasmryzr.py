#!/usr/bin/env python3

import sys,pandas as pd
import argparse, subprocess,pathlib
import pandas as pd
import numpy as np
from ast import literal_eval
import json 

def check_asm(contigs:list)-> tuple:

    q1, q3 = np.percentile(contigs, [25, 75])
    iqr = q3 - q1
    lower_bound = q1 - (1.5 * iqr)
    upper_bound = q3 + (1.5 * iqr)
    # print(f"Lower bound: {lower_bound}, Upper bound: {upper_bound}")
    return lower_bound, upper_bound

def check_species(species:list) -> str:
    """
    Check if the species is the same in the assembly and reads
    """
    sp = set(species)
    if len(sp) == 1:
        return 1
    elif len(sp) > 1:
        return f"Species from reads: {species[0]} and assembly: {species[1]} are different."
def check_val(val:float, min_val:float,metric) -> str:
    """
    Check if the coverage is above the minimum coverage
    """
    if val[1]:
        return 1
    if val[0] == "":
        return 1
    elif int(val[0]) >= min_val:
        return 1
    else:
        return f"{metric}: {val}, should be at least {min_val}`"

def check_contigs(upper, lower, contigs:list) -> str:
    """
    Check if the number of contigs is within the expected range
    """
    # print(contigs[0])
    # print(lower, upper)
    if contigs[1]:
        return 1
    if (contigs[0] >= lower) and (contigs[0] <= upper):
        return 1
    else:
        return f"Number of contigs: {contigs[0]} is outside the expected range of {lower} and {upper}."

def check_filesize(filesize) -> str:
    """
    Check if the file size is above the minimum size
    """
    if filesize[1] == "":
        return 1
    if filesize[0] == ">20000000":
        return 1
    else:
        return "File below recommended size of 20Mbp"

def make_comment(comments:list) -> str:
    """
    Make a comment from the list of comments
    """
    cmt = []
    for i in comments:
        if i != 1:
            cmt.append(i)
    if cmt:
        return " | ".join(cmt)
    return ""



def _generate_summary_table(results_files: list, output:list, min_depth:40, minquality : 30, minaln:70) -> list:
    print("Generating summary table")
    list_of_filename = {
        "read_assessment.txt" : ["Isolate","Reads","GC content", "Est average depth", "is_control", "filesize", "Qscore"],
        "assembly_assesment.txt":["Isolate","bp","# Contigs","N50"],
        "core_genome_stats.txt":["Isolate","% Aligned"],
        "speciation.txt":["Isolate","Species (reads)","Match 1 (reads)", "Match 1 (asm)"],
        "mlst.txt":["Isolate","Scheme","ST"],
        # "cluster.txt":[]
    }
    # print(list_of_filename)
    summary = pd.DataFrame()
    for file in results_files:
        # print(file)
        if pathlib.Path(file).name in list_of_filename and pathlib.Path(file).exists():
            df = pd.read_csv(file, sep = '\t')
            
            cols_to_keep = []
            for col in list_of_filename[pathlib.Path(file).name]:
                if col in df.columns:
                    cols_to_keep.append(col)
            df = df[cols_to_keep] 
            
            if summary.empty:
                summary = df
            else:
                summary = pd.merge(summary, df, how = 'outer', on = "Isolate")
    summary["is_control"] = summary["is_control"].fillna(False)
    summary = summary.fillna("")
    summary = summary.rename(columns = {"Match 1 (reads)":"Species (reads)","Match 1 (asm)":"Species (assembly)"})
    # summary["Data assessment"] = "ok"
    # summary["Comment"] = ""
    
    sp_cols = [i for i in summary.columns if "Species" in i]
    # print(summary)
    if sp_cols:
        summary["Species check"] = summary[sp_cols].apply(lambda x: check_species(x.tolist()), axis=1)
    if "filesize" in summary.columns:
        summary["File size check"] = summary[["filesize","Reads"]].apply(lambda x: check_filesize(x), axis=1)
        print(summary)
    if "Est average depth" in summary.columns:
        summary["Coverage check"] = summary[["Est average depth","is_control"]].apply(lambda x: check_val(x, min_depth, "Avg depth"),axis=1)
    if "% Aligned" in summary.columns:
        summary["Alignment check"] = summary[["% Aligned", "is_control"]].apply(lambda x:check_val(x, minaln, "Alignment %"), axis=1)

    if "# Contigs" in summary.columns:
        print("Checking number of contigs")
        bounds = check_asm(summary["# Contigs"].tolist())
        summary["Contigs check"] = summary[["# Contigs","is_control"]].apply(lambda x: check_contigs(bounds[1], bounds[0], x), axis=1)

    check_cols = [i for i in summary.columns if "check" in i]
    summary = summary.fillna("")
    summary["Comment"] = summary[check_cols].apply(lambda x: make_comment(x.tolist()), axis=1)
    summary["Data assessment"] = summary[check_cols].apply(lambda x: 1 if all(i == 1 for i in x.tolist()) else 0, axis=1)
    
    # print(summary)
    summary.drop(columns=check_cols, inplace=True)
    print(summary)
    if not summary.empty:
        summary.to_csv("summary.tsv", sep = '\t', index = False)
        output.append("summary.tsv")
        return "-f summary.tsv", output
    return "",output 


def _extract_tree(results_files: list, output : list) -> str:
    """
    Extract the tree from the results files
    """
    # check if the file exists
    
    for file in results_files:
        if pathlib.Path(file).exists() and "newick" in file:
            output.append(file)
            return f"-tr {file}",output
    return "",output

def _extract_distance_matrix(results_files: list, output:list) -> str:
    """
    Extract the distance matrix from the results files
    """
    # print(output)
    for file in results_files:
        if pathlib.Path(file).exists() and "distances" in file and "tsv" in file:       
            output.append(file)
            return f"-dm {file}", output
    # print(output)
    return "",output

def _extract_core_genome(results_files: list, output:list) -> str:
    """
    Extract the core genome from the results files
    """
    for file in results_files:
        if pathlib.Path(file).exists() and "vcf" in file:
            output.append(file)
            return f"-cg {file}",output
    return "",output

def _extract_core_genome_report(results_files: list, output:list) -> str:
    """
    Extract the core genome report from the results files
    """
    for file in results_files:
        if pathlib.Path(file).exists() and  "core_genome_stats" in file:
            output.append(file)
            return f"-cgr {file}",output
    return "",output

def _make_annotation_file(input_file: list, result_files:list, annot_cols : str) -> str:
    """

    Extract the annotation file from the input files

    """
    df = pd.DataFrame()
    annot_cols = annot_cols.split(",") if annot_cols else []
    # print(annot_cols)
    # print(pathlib.Path(input_file).exists())
    if len(annot_cols) > 0 and pathlib.Path(input_file).exists():
        # check if the file exists
        cols = ["Isolate"]
        # cols.extend(annot_cols)
        df = pd.read_csv(input_file, sep = '\t')
        for a in annot_cols:
            if a in df.columns:
                cols.append(a)
        df = df[cols]
        # print(df)
        
    for _file in result_files:
        # print(_file)
        if "cluster" in _file:
            # print("Found cluster file")
            clst = pd.read_csv(_file, sep = '\t')
            clst = clst.rename(columns = {"ID":"Isolate"})
            # print(clst)
            if not df.empty:
                df = df.merge(clst, how = 'left', on ="Isolate")
            else:
                df = clst
        else:
            try:
                tmp = pd.read_csv(_file, sep = '\t')
                # print(tmp)
                tmp_cols = ["Isolate"]
                for col in tmp.columns:
                    if col in annot_cols:
                        tmp_cols.append(col)
                # print(tmp_cols)
                if tmp_cols != ["Isolate"]:
                    tmp = tmp[tmp_cols]
                    # print(tmp)
                    # print(df)
                    if not df.empty:
                        df = df.merge(tmp, how = 'left', on ="Isolate")
                    else:
                        df = tmp
            except Exception as e:
                print(f"Error reading file {_file}: {e}")
        # print(df)
    
    if not df.empty:
        df.to_csv("annotation_file.tsv", sep = '\t', index = False)
    
    return   "--annotate annotation_file.tsv" if not df.empty else ""
    
def _get_reference(reference: str) -> str:
    """
    Extract the reference from the results files
    """
    if reference != "" and pathlib.Path(reference).exists():
        return f"-r {reference}"
    return ""

def _get_mask(mask: str) -> str:
    """
    Extract the mask from the results files
    """
    if mask != "" and pathlib.Path(mask).exists():
        return f"-m {mask}"
    return ""

def _get_other_files(results_files: list,ouput:list) -> str:
    """
    Extract the other files from the results files
    """
    other_files = []
    for file in results_files:
        # print(file)
        if file not in ouput and pathlib.Path(file).exists() and ("txt" in file or "tsv" in file or "json" in file):
            other_files.append(f"-f {file}")
        

    return " ".join(other_files) if other_files else ""

def _extract_pangenome_classification(results_files: list, output:list) -> str:
    """
    Extract the pangenome summary from the results files
    """
    for file in results_files:
        if pathlib.Path(file).exists() and "classification.tab" in file:
            output.append(file)
            return f"--pangenome_characterization {file}",output
    return "",output

def _extract_pangenome_rtab(results_files: list, output:list) -> str:
    """
    Extract the pangenome rtab from the results files
    """
    for file in results_files:
        if pathlib.Path(file).exists() and "gene_presence_absence.Rtab" in file:
            output.append(file)
            return f"--pangenome_rtab {file}",output
    return "",output

def _extract_pangenome_groups(results_files: list, output:list) -> str:
    """
    Extract the pangenome groups from the results files
    """
    for file in results_files:
        # print(file)
        if pathlib.Path(file).exists() and "group.txt" in file:
            output.append(file)
            return f"--pangenome_groups {file}",output
    return "",output

def _combine_reads_iqr(results_files: list, output:list) -> str:

    """
    Combine the reads iqr from the results files
    """
    read_assessment = pd.DataFrame()
    for file in results_files:
        print(file)
        if pathlib.Path(file).exists() and pathlib.Path(file).name in ["read_assessment.txt","insertiqr.txt"]:
            iiqr = pd.read_csv(file, sep = '\t')
            print(iiqr)
            if read_assessment.empty:
                read_assessment = iiqr
            else:
                read_assessment = pd.merge(read_assessment, iiqr, how = 'outer', on = "Isolate")
            print(f"Should be removing file: {file}")
            # results_files.remove(file)
   
        # if pathlib.Path(file).exists() and "read_assessment.txt" in file:
        #     ra = pd.read_csv(file, sep = '\t')
        #     print(ra)
        #     if read_assessment.empty:
        #         read_assessment = ra
        #     else:
        #         read_assessment = pd.merge(read_assessment, ra, how = 'outer', on = "Isolate")
        #     print(f"Should be removing file: {file}")
        #     results_files.remove(file)
    # print(read_assessment)
    if not read_assessment.empty:
            for f in ["read_assessment.txt", "insertiqr.txt"]:
                for i in results_files:
                    if f in i:
                        print(f"Removing {i} from results_files")
                        results_files.remove(i)
            
            read_assessment.to_csv("read_assessment.txt", sep = '\t', index = False)
            results_files.append("read_assessment.txt")
            output.append("read_assessment.txt")
            # print(results_files)
            return f"-f read_assessment.txt",output,results_files
    return "",output,results_files

def _run_datasmryzr(tree:str,
                    distance_matrix:str,
                    core_genome:str,
                    core_genome_report:str,
                    other_files:str,
                    reference,mask:str, 
                    annotation:str, 
                    bkgd_color:str,
                    text_color:str,
                    job_id:str,
                    pangenome_classification:str,
                    pangenome_rtab:str,
                    pangenome_groups:str,
                    read_assessment) -> str:    
    """
    Run the datasmryzr pipeline
    """
    cmd = f"datasmryzr -c bohra_config.json -bg '{bkgd_color}' -fc '{text_color}' {other_files} {pangenome_classification} {pangenome_rtab} {pangenome_groups} {tree} {distance_matrix} {core_genome} {core_genome_report} {reference} {mask} {annotation}"
    # print(cmd)
    # p = subprocess.run(cmd, shell=True, capture_output=True)
    # if p.returncode != 0:
    #     print(p.stderr.decode())
    #     raise Exception(f"Error running datasmryzr: {p.stderr.decode()}")
    # else:
    #     print(p.stdout.decode())
    #     print("datasmryzr run complete")
    #     return p.stdout.decode()

def generate_config(cluster_method:str,
                    cluster_threshold:str,
                    pangenome_groups:str,
                    kraken2_db:str,
                    speciation: str) -> None:
    
    key = {
        "GRP":pangenome_groups,
        "THRESHOLDS":cluster_threshold,
        "METHOD":cluster_method,
        "KRAKEN2_DB": kraken2_db,}
    
    with open(f"{pathlib.Path(__file__).parent / 'base_config.json'}","r") as j:
        cfg = j.read()
        
        for rpl in key:
            if key[rpl] != "":
                cfg = cfg.replace(rpl, key[rpl])
        # print(cfg)
        dct = literal_eval(cfg)
        if "speciation" in dct["comments"]:
            sp = "kraken2" if speciation == "kraken2" else "sylph"
            sprow = dct["comments"]["speciation"].get(sp, [])
            dct["comments"]["speciation"] = sprow
        print(dct["comments"]["speciation"])
    with open("bohra_config.json" , "w") as o:
        json.dump(dct, o , indent = 4)

        

    


def _compile(args):
    # print(f"{args.annot_cols}")
    # print(args.results_files)
    # print(args.speciation)
    output = []
    results_files = args.results_files
    results_files = [i for i in results_files if pathlib.Path(i).exists()]
    read_assessment,output,results_files = _combine_reads_iqr(results_files, output)
    tree,output = _extract_tree(results_files, args.job_id)
    distance_matrix,output = _extract_distance_matrix(results_files, output)
    core_genome,output = _extract_core_genome(results_files, output)
    core_genome_report,output = _extract_core_genome_report(results_files, output)
    summary,output = _generate_summary_table(results_files, output, 40, 30, 70)
    print(output)
    
    panclass,output = _extract_pangenome_classification(results_files, output)
    panrtab,output= _extract_pangenome_rtab(results_files, output)
    pangroups,output = _extract_pangenome_groups(results_files, output)
    # print(panclass, panrtab, pangroups)
    other_files = _get_other_files(results_files,output)
    other_files = other_files + " " + summary if summary else other_files
    reference = _get_reference(args.reference)
    mask = _get_mask(args.mask)
    annotation = _make_annotation_file(args.input_file, results_files, f"{args.annot_cols}")
    generate_config(args.cluster_method, args.cluster_threshold, args.pangenome_groups, args.kraken2_db, 'kraken2' if args.speciation == 'true' else "sylph")
    p = _run_datasmryzr(tree,
                        distance_matrix,
                        core_genome,
                        core_genome_report,
                        other_files,
                        reference,
                        mask,
                        annotation,
                        args.bkgd,
                        args.text_color,
                        args.job_id,
                        panclass,
                        panrtab,
                        pangroups,
                        read_assessment)
    


def set_parsers():
    # setup the parser
    parser = argparse.ArgumentParser(description='Collate and write bohra report',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--mask',
        help='',
        default = '')
    parser.add_argument('--job_id',
        help='',
        default = '',
        nargs='+')
    parser.add_argument('--results_files',
        help = '',
        default = '',
        nargs='+')
    parser.add_argument('--reference',
        help = '',
        default = '')
    parser.add_argument('--input_file',
        help = '',
        default = '')
    parser.add_argument("--annot_cols",
        help="",
        default = '',
    )
    parser.add_argument('--bkgd',
        help='',
        default = '')
    parser.add_argument('--text_color',
        help='',
        default = '',)
    parser.add_argument('--cluster_method',
                        help = '',
                        default = ''
                        )
    parser.add_argument('--cluster_threshold',
    help = '',
    default = ''
    )
    parser.add_argument('--pangenome_groups',
    help = '',
    default = ''
    )
    parser.add_argument('--speciation',
    help = '',
    default = ''
    )
    parser.add_argument('--kraken2_db',
    help = '',
    default = ''
    )
 
    
    parser.set_defaults(func=_compile)
    args = parser.parse_args()
    
    if vars(args) == {}:
        parser.print_help(sys.stderr)
    else:
        # print(args)
        args.func(args)
	

def main():
    """
    run pipeline
    """

    args = set_parsers()
    

if __name__ == "__main__":
    main()

