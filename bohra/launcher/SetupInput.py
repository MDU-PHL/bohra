import pathlib
import datetime
import pandas as pd
import logging
from bohra.launcher.Utils import CustomFormatter, _check_path
# Logger
# Logger
LOGGER =logging.getLogger(__name__) 
LOGGER.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(CustomFormatter())
fh = logging.FileHandler('bohra.log')
fh.setLevel(logging.DEBUG)
formatter = logging.Formatter('[%(levelname)s:%(asctime)s] %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p') 
fh.setFormatter(formatter)
LOGGER.addHandler(ch) 
LOGGER.addHandler(fh)


CFG = {
    "reads": {
        "ext":"*.f*q.gz",
        "num_expected":2,
    },
    "se-reads": {
        "ext":"*.f*q*",
        "num_expected":1
    },
    "ont": {
        "ext":"*.f*q*",
        "num_expected":1
    },
    "asm": {
        "ext":"*.f*a",
        "num_expected":1
    }
}
  
    
def _get_isolate_list(reads:str,
                      contigs:str,
                      ) -> list:
    """
    Retrieve a list of unique isolate names from files in the specified directory.

    :param _dir: Path to the directory containing the files.
    :param isolates: A list of isolate names to filter the results. If empty, all isolates are returned.
    :param ext: File extension pattern to search for (default is "*.f*q.gz").
    :return: A list of unique isolate names found in the directory.
    """
    # iso_found = _get_isolate_list(_dir=_dir, isolates=isolates, ext="*.f*q.gz")
    all_data = []
    if reads != '' and _check_path(reads):
        ext = CFG['reads']['ext']
        reads_list = sorted(pathlib.Path(reads).rglob(ext))
        all_data.extend(reads_list)
    elif contigs != '' and _check_path(contigs):
        ext = CFG['asm']['ext']
        contigs_list = sorted(pathlib.Path(contigs).rglob(ext))
        all_data.extend(contigs_list)
    else:
        LOGGER.error("No valid reads or contigs path provided. Please check your input.")
        raise SystemExit
    
    iso_found = set()

    for sq in all_data:
        nme = sq.name.split('_')[0]
        iso_found.add(nme)

    return list(iso_found)

def _order_pereads(reads: list) -> list:
    """
    Order paired-end reads by their names.

    :param reads: List of read file paths.
    :return: Ordered list of read file paths.
    """
    
    
    r1 = [r for r in reads if "_R1" in r.name or "_1.f" in r.name or "_r1" in r.name]
    r1 = sorted(r1, key=lambda f: f.stat().st_ctime)[-1]
    r2 = [r for r in reads if "_R2" in r.name or "_2.f" in r.name or "_r2" in r.name]
    r2 = sorted(r2, key=lambda f: f.stat().st_ctime)[-1]
    
    return [r1,r2]

def _glob_sequences(_dir: pathlib.Path, 
                    isolates: list, 
                    sequence_type: str) -> pd.DataFrame:
    """
    Process and validate sequence files for the specified sequence type.

    :param _dir: Path to the directory containing sequence files.
    :param isolates: A list of isolate names to filter the results.
    :param sequence_type: Type of sequence data to process (e.g., 'reads', 'asm').
    :return: True if processing is successful, otherwise raises an error.
    """
    data = []
    if sequence_type in CFG:
        
        for iso in isolates:
            
            seqs = sorted(pathlib.Path(_dir).rglob(f"*{iso}{CFG[sequence_type]['ext']}"))
            LOGGER.info(f"Found {len(seqs)} {sequence_type} for {iso} in {_dir}.")
            if seqs == []:
                seqs = sorted(pathlib.Path(_dir).rglob(f"{iso}/{CFG[sequence_type]['ext']}"))
            if len(seqs) >= CFG[sequence_type]['num_expected']:
                LOGGER.info(f"Now add {sequence_type} for {iso}")
                if sequence_type == "reads":
                    seqs = _order_pereads(seqs)
                    heads = [f"r{i+1}" for i in range(len(seqs))]
                    if len(seqs) > CFG[sequence_type]['num_expected']:
                        LOGGER.warning(f"Expected 2 reads for {iso}, but found {len(seqs)}. Will use the most recent.")
                elif sequence_type == "asm":
                    heads = ["assembly"]
                else:
                    LOGGER.error(f"Unexpected sequence type: {sequence_type}.")
                    raise SystemExit
                line = dict(zip(heads,seqs))
                line['Isolate'] = iso
                data.append(line)
            elif len(seqs) < CFG[sequence_type]['num_expected']:
                LOGGER.warning(f"There do not appear to be {CFG[sequence_type]['num_expected']} {sequence_type} for {iso}. Skipping")
            
        if data != []:
            return pd.DataFrame(data)
        else:
            LOGGER.warning(f"It appears that no {sequence_type} have been found. Please check your input and try again.")
            raise SystemExit
        
    else:
        LOGGER.critical(f"An error has occurred somewhere - {sequence_type} is not a valid type of input data.")
        raise SystemExit

        
def _extract_isolates(isolate_ids:str)-> pd.DataFrame:
    """
    Extract isolates from isolate_ids file
    :input - path to isolate_ids file
    :output - list of isolates
    """
    try:
        isolates_df = pd.read_csv(isolate_ids, sep="\t")
        if 'Isolate' not in isolates_df.columns:
            LOGGER.error("The provided isolate_ids file does not contain a column named 'Isolate'. Please check the file format.")
            raise SystemExit
        if 'species' not in isolates_df.columns:
            LOGGER.warning("The provided isolate_ids file does not contain a column named 'species'. This column is optional but recommended for further annotation of trees.")
        
    except Exception as e:
        LOGGER.error(f"An error occurred while reading the isolate_ids file: {e}")
        raise SystemExit
    
    return isolates_df
    
def find_data(reads:str,
              contigs:str,
              isolate_ids:str,
              outname:str,
              ) -> None:
    """
    Find data in the path provided
    :input - path to reads/contigs
                isolate_ids - file containing isolate ids
                input_type - type of input data (reads, asm)
    :output - list of data
    """
    isolates = []
    bohra_table = []
    if isolate_ids != '' and _check_path(isolate_ids):
        LOGGER.info(f"Extracting isolates from {isolate_ids}")
        isolates_df = _extract_isolates(isolate_ids)
        bohra_table.append(isolates_df)
        isolates = list(isolates_df['Isolate'].unique())
    else:
        LOGGER.warning("No isolate_ids file provided. Using all isolates found in the reads/contigs directories.")
        isolates = _get_isolate_list(reads=reads, contigs = contigs)
    # now get reads
    if reads != '' and _check_path(reads):
        reads_path = pathlib.Path(reads)
        if not reads_path.exists():
            LOGGER.error(f"The specified reads path {reads} does not exist.")
            raise SystemExit
        reads_tab = _glob_sequences(_dir=reads_path, isolates=isolates, sequence_type='reads')
        if not reads_tab.empty:
            bohra_table.append(reads_tab)
    if contigs != '' and _check_path(contigs):
        contigs_path = pathlib.Path(contigs)
        if not contigs_path.exists():
            LOGGER.error(f"The specified contigs path {contigs} does not exist.")
            raise SystemExit
        contigs_tab = _glob_sequences(_dir=contigs_path, isolates=isolates, sequence_type='asm')
        if not contigs_tab.empty:
            bohra_table.append(contigs_tab)
    
    if bohra_table:

        bohra_df = pd.DataFrame()
        for df in bohra_table:
            if bohra_df.empty:
                bohra_df = df
            else:
                bohra_df = pd.merge(bohra_df, df, on='Isolate', how='outer')
        # bohra_df = bohra_df.loc[:,~bohra_df.columns.duplicated()]
        bohra_df = bohra_df.sort_values(by='Isolate')
        bohra_df = bohra_df.fillna('not_supplied')
        bohra_df.to_csv(outname,sep = "\t", index=False)
        LOGGER.info(f"Input files generated successfully and saved to '{outname}'.")