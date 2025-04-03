import pathlib
import datetime
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
    "pe-reads": {
        "ext":"*.f*q.gz",
        "num_expected":2
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
  
    
def _get_isolate_list(_dir: pathlib.Path,
                      isolates: list, 
                      ext: str = f"*.f*q.gz") -> list:
    """
    Retrieve a list of unique isolate names from files in the specified directory.

    :param _dir: Path to the directory containing the files.
    :param isolates: A list of isolate names to filter the results. If empty, all isolates are returned.
    :param ext: File extension pattern to search for (default is "*.f*q.gz").
    :return: A list of unique isolate names found in the directory.
    """
    iso_found = _get_isolate_list(_dir=_dir, ext="*.f*q.gz")
    
    all_data = sorted(_dir.rglob(ext))

    iso_found = set()

    for reads in all_data:
        nme = reads.name.split('_')[0]

        iso_found.add(nme)

    isolist = []
    if isolates != []:

        for i in list(iso_found):
            for j in isolates:
                if i in j:
                    isolist.append(i)
    else:
        isolist = iso_found
    
    LOGGER.info(f"Found {len(list(iso_found))} distinct sample names at path {_dir}")
    return list(iso_found)

def _glob_sequences(_dir: pathlib.Path, 
                    isolates: list, 
                    sequence_type: str) -> bool:
    """
    Process and validate sequence files for the specified sequence type.

    :param _dir: Path to the directory containing sequence files.
    :param isolates: A list of isolate names to filter the results.
    :param sequence_type: Type of sequence data to process (e.g., 'reads', 'asm').
    :return: True if processing is successful, otherwise raises an error.
    """
    if sequence_type in CFG:
        isolist = _get_isolate_list(_dir=_dir, isolates=isolates, ext=CFG[sequence_type]["ext"])

        lines = []

        for iso in isolist:
            reads = sorted(_dir.rglob(f"*{iso}*.f*q.gz"))
            if len(reads) == CFG[sequence_type]['num_expected']:
                LOGGER.info(f"Now add reads for {iso}")
                lines.append(f"{iso}\t{reads[0]}\t{reads[1]}")
            elif len(reads) < CFG[sequence_type]['num_expected']:
                LOGGER.warning(f"There do not appear to be {CFG[sequence_type]['num_expected']} reads for {iso}. Skipping")
            else:
                LOGGER.warning(f"There appear to be more than {CFG[sequence_type]['num_expected']} reads available for {iso}. Skipping")
        if lines != []:
            LOGGER.info(f"Saving reads file as isolates.tab")
            pathlib.Path(f"{sequence_type}.tab").write_text('\n'.join(lines))
        else:
            LOGGER.warning(f"It appears that no reads have been found. Please check your input and try again.")
        
        return True
    else:
        LOGGER.critical(f"An error has occurred somewhere - {sequence_type} is not a valid type of input data.")
        raise SystemExit

        
def _extract_isolates(isolate_ids:str):
    """
    Extract isolates from isolate_ids file
    :input - path to isolate_ids file
    :output - list of isolates
    """
    with open(isolate_ids, 'r') as f:
        isolates = f.read().strip().split('\n')
    
    return isolates

    
def find_data( reads:bool,
                contigs:bool,
                isolate_ids:str,
                path:str):
    """
    Find data in the path provided
    :input - path to reads/contigs
                isolate_ids - file containing isolate ids
                reads - boolean
                contigs - boolean
    :output - list of data
    """
    LOGGER.info(f"Finding data in {path}")
    if path != '' and _check_path(path):
        LOGGER.info(f"Path provided : {path} exists.")
        if isolate_ids != '': 
            isolates = _extract_isolates(isolate_file= isolate_ids) if _check_path(isolate_ids) else []
        else:
            isolates = []
        _glob_sequences(path = path, isolates= isolates, sequence_type = 'reads' if reads else 'asm')
    else:
        LOGGER.critical(f"There seems to be a problem with your read path. Please provide a valid path to the reads you wish to include")
        raise SystemExit
