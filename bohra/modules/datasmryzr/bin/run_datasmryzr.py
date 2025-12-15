#!/usr/bin/env python3

import sys,pandas as pd
import argparse, subprocess,pathlib
import pandas as pd
import numpy as np
from ast import literal_eval
import json 

def check_asm(contigs:list)-> tuple:
    contigs = [i for i in contigs if not isinstance(i, str)]
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
    sp = set([i for i in species if i != ""])
    if len(sp) == 1:
        return 1
    elif len(sp) > 1:
        return f"Species inconsistency."
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
        return f"{metric}: {val[0]}, should be at least {min_val}`"

def check_val_aln(val:float, min_val:float,metric:str,strict:bool) -> str:
    """
    Check if the coverage is above the minimum coverage and not an outlier
    """
    if val[1]:
        return 1
    if val[0] == "" and val[2] == "":
        return 1
    elif int(val[0]) >= min_val and val[2] == "":
        return 1
    elif val[2] != "":
        return f"Isolate marked as outlier in core genome analysis and is remove from comparative dataset." if strict else f"Isolate marked as outlier in core genome analysis. Should be removed from comparative dataset and bohra rerun." 
    else:
        return f"{metric}: {val[0]}, should be at least {min_val}`"



def check_contigs(upper, lower, contigs:list) -> str:
    """
    Check if the number of contigs is within the expected range
    """
    # print(contigs[0])
    # print(lower, upper)
    if contigs[1]:
        return 1
    if contigs[0] == "":
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



def _generate_summary_table(input_file: str, results_files: list, output:list, min_depth:40, minquality : 30, minaln:70, strict:str) -> list:
    print("Generating summary table")
    list_of_filename = {
        "read_assessment.txt" : ["Isolate","Reads","GC", "Depth","Genome size", "is_control", "filesize", "Qscore"],
        "assembly_assessment.txt":["Isolate","Length","# Contigs","Assembly N50"],
        "core_genome_stats.txt":["Isolate","% Aligned", "Aln_outlier"],
        "speciation.txt":["Isolate","Species (reads)","Match 1 (reads)", "Match 1 (asm)"],
        "mlst.txt":["Isolate","Scheme","ST"],
        # "cluster.txt":[]
    }
    strict = True if strict.lower() == "false" else False
    colsblcklst = [ "r1","r2","assembly", "is_control"]
    input_data = pd.read_csv(input_file, sep='\t')
    input_data = input_data.rename(columns = {"species":"Species_expected"})
    input_data = input_data.replace("not_supplied","")
    for col in input_data.columns:
        if col in colsblcklst:
            input_data = input_data.drop(columns=[col])
    # print(list_of_filename)
    summary = pd.DataFrame()
    for file in results_files:
        # print(file)
        if pathlib.Path(file).name in list_of_filename and pathlib.Path(file).exists():
            if 'mlst' in pathlib.Path(file).name:
                df = pd.read_csv(file, sep = '\t', dtype=str)
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
    summary["Comment"] = ""
    summary["Data assessment"] = 1
    summary["is_control"] = summary["is_control"].fillna(False)
    summary = pd.merge(input_data, summary, how = 'left', on = "Isolate")
    summary = summary.fillna("")
    summary = summary.rename(columns = {"Match 1 (reads)":"Species (reads)","Match 1 (asm)":"Species (assembly)"})
    if 'ST' in summary.columns:
        summary['ST'] = summary['ST'].apply(lambda x: str(int(x)) if not isinstance(x, str) else x)
    int_cols = ["Reads","Assembly length","# Contigs","Assembly N50", "Depth", "Genome size"]
    for i in int_cols:
        if i in summary.columns:
            summary[i] = summary[i].apply(lambda x: int(x) if x != "" else "")
    
    sp_cols = [i for i in summary.columns if "Species" in i]
    # print(summary)
    if sp_cols:
        summary["Species check"] = summary[sp_cols].apply(lambda x: check_species(x.tolist()), axis=1)
        summary["Species"] = summary[sp_cols].apply(lambda x: ':'.join(list(set(x))).strip(":"), axis=1)
        summary.drop(columns=sp_cols, inplace=True)
        cols = [i for i in summary.columns if i not in sp_cols]
    if "filesize" in summary.columns:
        summary["File size check"] = summary[["filesize","Reads"]].apply(lambda x: check_filesize(x), axis=1)
        print(summary)
    if "Depth" in summary.columns:
        summary["Coverage check"] = summary[["Depth","is_control"]].apply(lambda x: check_val(x, min_depth, "Avg depth"),axis=1)
    if "% Aligned" in summary.columns:
        summary["Alignment check"] = summary[["% Aligned", "is_control", "Aln_outlier"]].apply(lambda x:check_val_aln(x, minaln, "Alignment %",strict), axis=1)

    if "# Contigs" in summary.columns:
        print("Checking number of contigs")
        bounds = check_asm(summary["# Contigs"].tolist())
        summary["Contigs check"] = summary[["# Contigs","is_control"]].apply(lambda x: check_contigs(bounds[1], bounds[0], x), axis=1)

    check_cols = [i for i in list(summary.columns) if "check" in i]
    print(check_cols)
    check_cols.append("File size check")
    
    print(summary.columns)
    summary = summary.fillna("")
    print(summary)
    summary["Comment"] = summary[check_cols].apply(lambda x: make_comment(x), axis=1)
    summary["Data assessment"] = summary[check_cols].apply(lambda x: 1 if all(i == 1 for i in x.tolist()) else 0, axis=1)
    if "Aln_outlier" in summary.columns:
        check_cols.append("Aln_outlier")
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

def _extract_cluster_table(results_files: list, output:list) -> str:
    """
    Extract the distance matrix from the results files
    """
    # print(output)
    for file in results_files:
        if pathlib.Path(file).exists() and "clusters" in file:
            try:
                pd.read_csv(file, sep = '\t')       
                output.append(file)
                return f"-ct {file}", output
            except Exception as e:
                print(f"Error reading cluster file {file}: {e}")
                continue
    # print(output)
    return "",output

def _extract_distance_matrix(results_files: list, output:list) -> str:
    """
    Extract the distance matrix from the results files
    """
    # print(output)
    for file in results_files:
        if pathlib.Path(file).exists() and "distances" in file and "tsv" in file:       
            try:
                pd.read_csv(file, sep = '\t')       
                output.append(file)
                return f"-dm {file}", output
            except Exception as e:
                print(f"Error reading cluster file {file}: {e}")
                continue
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
                    cluster_table:str,
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
                    read_assessment,
                    pipeline:str,
                    pipeline_version:str,
                    no_downloadable_tables:str) -> str:    
    """
    Run the datasmryzr pipeline
    """
    cmd = f"datasmryzr --title '{job_id}' -c bohra_config.json -bg '{bkgd_color}' -fc '{text_color}' --pipeline {pipeline} --pipeline_version '{pipeline_version}' {no_downloadable_tables} {other_files} {pangenome_classification} {pangenome_rtab} {pangenome_groups} {tree} {distance_matrix} {cluster_table} {core_genome} {core_genome_report} {reference} {mask} {annotation} {read_assessment}"
    print(cmd)
    p = subprocess.run(cmd, shell=True, capture_output=True)
    if p.returncode != 0:
        print(p.stderr.decode())
        raise Exception(f"Error running datasmryzr: {p.stderr.decode()}")
    else:
        print(p.stdout.decode())
        print("datasmryzr run complete")
        return p.stdout.decode()

def extract_cmd(launchdir: str) -> str:
    """
    Extract the command from the launch directory
    """
    cmd = ""
    if pathlib.Path(launchdir).exists():
        with open(pathlib.Path(launchdir) / "bohra_run.log", "r") as f:
            cmd = f.read().strip().split("\n")[-1].strip('\x1b[1m[0m')
        if cmd.startswith("nextflow"):
            # cmd = cmd.replace("nextflow", "bohra")
            return cmd
        else:
            return ""
        
def get_ref_accession(reference: str) -> str:
    """
    Get the accession number from the reference file
    """
    if pathlib.Path(reference).exists():
        acc = ""
        with open(reference, "r") as f:
            header = f.readline().strip().split("\n")[0]

        if header.startswith(">"):
            acc = header.strip(">")
        elif header.startswith("LOCUS"):
            acc = header.split(" ")[1]
        return acc
    return ""

def generate_config(cluster_method:str,
                    cluster_threshold:str,
                    pangenome_groups:str,
                    kraken2_db:str,
                    speciation: str,
                    launchdir: str,
                    reference:str) -> None:

    key = {
        "GRP":pangenome_groups,
        "THRESHOLDS":cluster_threshold,
        "METHOD":cluster_method,
        "KRAKEN2_DB": kraken2_db,
        "CMD": extract_cmd(launchdir),
        "REF_FILE": reference if pathlib.Path(reference).exists() else "",
        "REF_ACCESSION": get_ref_accession(reference) if pathlib.Path(reference).exists() else "",
        }

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

        

    
def get_pipeline_version(results_files: list) -> str:
    """
    Get the pipeline version from the results files
    """
    for file in results_files:
        try:
            if pathlib.Path(file).exists() and "version" in file:
                versions = pd.read_csv(file, sep = '\t')
                version = versions[versions['tool'] == 'bohra']['version'].values
                return version[0] if len(version) > 0 else "not provided"
        except Exception as e:
            print(f"Error reading file {file}: {e}")
            return "not provided"

def _compile(args):
    # print(f"{args.annot_cols}")
    # print(args.results_files)
    # print(args.speciation)
    output = []
    results_files = args.results_files
    results_files = [i for i in results_files if pathlib.Path(i).exists()]
    read_assessment,output,results_files = _combine_reads_iqr(results_files, output)
    tree,output = _extract_tree(results_files, output)
    distance_matrix,output = _extract_distance_matrix(results_files, output)
    cluster_table,output = _extract_cluster_table(results_files, output)
    core_genome,output = _extract_core_genome(results_files, output)
    core_genome_report,output = _extract_core_genome_report(results_files, output)
    summary,output = _generate_summary_table(args.input_file, results_files, output, 40, 30, 70,args.ignore_warnings)
    print(read_assessment)
    
    panclass,output = _extract_pangenome_classification(results_files, output)
    panrtab,output= _extract_pangenome_rtab(results_files, output)
    pangroups,output = _extract_pangenome_groups(results_files, output)
    # print(panclass, panrtab, pangroups)
    other_files = _get_other_files(results_files,output)
    other_files = other_files + " " + summary if summary else other_files
    reference = _get_reference(args.reference)
    mask = _get_mask(args.mask)
    annotation = _make_annotation_file(args.input_file, results_files, f"{args.annot_cols}")
    generate_config(args.cluster_method, args.cluster_threshold, args.pangenome_groups, args.kraken2_db, 'kraken2' if args.speciation == 'true' else "sylph", args.launchdir, args.reference)
    pipeline_version = get_pipeline_version(args.results_files)
    ndt = '--no-downloadable-tables' if args.no_downloadable_tables.lower() == 'true' else ''
    p = _run_datasmryzr(tree,
                        distance_matrix,
                        cluster_table,
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
                        read_assessment,
                        args.pipeline,
                        pipeline_version,
                        ndt
                        )
    


def set_parsers():
    # setup the parser
    parser = argparse.ArgumentParser(description='Collate and write bohra report',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--mask',
        help='',
        default = '')
    parser.add_argument('--job_id',
        help='',
        default = '')
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
    parser.add_argument("--ignore_warnings",
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
    parser.add_argument('--launchdir',
    help = '',
    default = ''
    )
    parser.add_argument('--pipeline',
    help = '',
    default = ''
    )
    parser.add_argument('--report_outdir',
    help = '',
    default = ''
    )
    parser.add_argument('--no-downloadable-tables',
        help='Disable downloadable tables in the report html.',
        default='false'
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

