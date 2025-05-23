#!/usr/bin/env python3

import sys,pandas as pd
import argparse, subprocess,pathlib


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
    print(output)
    for file in results_files:
        if pathlib.Path(file).exists() and "distances" in file and "tsv" in file:       
            output.append(file)
            return f"-dm {file}", output
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

def _make_annotation_file(input_file: list, result_files:list, annot_cols : list) -> str:
    """

    Extract the annotation file from the input files

    """
    df = pd.DataFrame()
    print(input_file)
    if len(annot_cols) > 0 and pathlib.Path(input_file).exists():
        # check if the file exists
        cols = ["Isolate"]
        cols.extend(annot_cols)
        df = pd.read_csv(input_file, sep = '\t', usecols = cols)

        print(df)
    cluster_path = [i for i in result_files if "cluster" in i]
    if len(cluster_path) != []:
        cluster_path = cluster_path[0]
        clst = pd.read_csv(cluster_path, sep = '\t')
        clst = clst.rename(columns = {"ID":"Isolate"})
        print(clst)
        if not df.empty:
            df = df.merge(clst, how = 'left', on ="Isolate")
        else:
            df = clst
    if not df.empty:
        df.to_csv("annotation_file.tsv", sep = '\t', index = False)
    
    return   "-a annotation_file.tsv" if not df.empty else ""
    
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
        if file not in ouput and pathlib.Path(file).exists() and ("txt" in file or "tsv" in file or "json" in file):
            other_files.append(f"-f {file}")
        

    return " ".join(other_files) if other_files else ""

def _run_datasmryzr(tree:str,
                    distance_matrix:str,
                    core_genome:str,
                    core_genome_report:str,
                    other_files:str,
                    reference,mask:str, 
                    annotation:str, 
                    bkgd_color:str,
                    text_color:str) -> str:    
    """
    Run the datasmryzr pipeline
    """
    cmd = f"datasmryzr -bg '{bkgd_color}' -fc '{text_color}' {other_files} {tree} {distance_matrix} {core_genome} {core_genome_report} {reference} {mask} {annotation}"
    print(cmd)
    p = subprocess.run(cmd, shell=True, capture_output=True)
    if p.returncode != 0:
        print(p.stderr.decode())
        raise Exception(f"Error running datasmryzr: {p.stderr.decode()}")
    else:
        print(p.stdout.decode())
        print("datasmryzr run complete")
        return p.stdout.decode()

def _compile(args):

    output = []
    tree,output = _extract_tree(args.results_files, output)
    distance_matrix,output = _extract_distance_matrix(args.results_files, output)
    core_genome,output = _extract_core_genome(args.results_files, output)
    core_genome_report,output = _extract_core_genome_report(args.results_files, output)
    other_files = _get_other_files(args.results_files,output)
    reference = _get_reference(args.reference)
    mask = _get_mask(args.mask)
    annotation = _make_annotation_file(args.input_file, args.results_files, args.annot_cols.split(","))

    p = _run_datasmryzr(tree,
                        distance_matrix,
                        core_genome,
                        core_genome_report,
                        other_files,
                        reference,
                        mask,
                        annotation,
                        args.bkgd,
                        args.text_color)
    


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
        default = '')
    
    
    
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

