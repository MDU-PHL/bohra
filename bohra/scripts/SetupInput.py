import pathlib
import datetime
import logging
from bohra.scripts.Utils import CustomFormatter, _check_path
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


class SetupInputFiles(object):
    def __init__(self, 
                 
                 ):

        self.read_path = args.read_path
        self.isolate_list = args.isolate_list
        self.now = datetime.datetime.today().strftime("%d_%m_%y_%H")
        self.day = datetime.datetime.today().strftime('%Y-%m-%d')

    def _check_path(self,path):

        """
        Check if path provided exists
        :input - path to where reads/contigs/inputfile are
        :output - boolean 
        """
        LOGGER.info(f"Checking for {path}")
        if pathlib.Path(path).exists():
            LOGGER.info(f"{path} exists.")
            return True
        else:
            LOGGER.critical(f"Path provided : {path} does not exist - please try again.")
            return False
        
   
    
    def _get_isolate_list(self, _dir,  ext= f"*.f*q.gz"):

        all_data = sorted(_dir.rglob(ext))

        iso_found = set()

        for reads in all_data:
            nme = reads.name.split('_')[0]

            iso_found.add(nme)
        
        LOGGER.info(f"Found {len(list(iso_found))} distinct sample names at path {_dir}")
        return list(iso_found)
    
    

def _glob_reads(self, _dir, isolates):

    iso_found = self._get_isolate_list(_dir = _dir, ext = "*.f*q.gz")
    
    isolist = []
    if isolates != []:

        for i in list(iso_found):
            for j in isolates:
                if i in j:
                    isolist.append(i)
    else:
        isolist = iso_found

    lines = []

    for iso in isolist:

        reads = sorted(_dir.rglob(f"*{iso}*.f*q.gz"))
        if len(reads) == 2:
            LOGGER.info(f"Now add reads for {iso}")
            lines.append(f"{iso}\t{reads[0]}\t{reads[1]}")
        elif len(reads) <2:
            LOGGER.warning(f"There do not appear to be 2 reads for {iso}. Skipping")
        else:
            LOGGER.warning(f"There appear to be more than 2 reads available for {iso}. Skipping")
    if lines != []:
        LOGGER.info(f"Saving reads file as isolates.tab")
        pathlib.Path('isolates.tab').write_text('\n'.join(lines))
    else:
        LOGGER.warning(f"It appears that no reads have been found. Please check your input and try again.")
    
    return True

        
def _glob_data(path, isolates =[], data_type = 'reads'):

    """
    glob path provided for data
    :input - path provided.
                data type (reads or contigs - default to reads)
                
    :output - list of data
    """

    _dir = pathlib.Path(path).resolve()

    if data_type == 'reads':

        _glob_reads(_dir = _dir, isolates= isolates)

def _extract_isolates(isolate_ids):
    """
    Extract isolates from isolate_ids file
    :input - path to isolate_ids file
    :output - list of isolates
    """
    with open(isolate_ids, 'r') as f:
        isolates = f.read().strip().split('\n')
    
    return isolates

    
def find_data( reads,
                contigs,
                isolate_ids,
                path):


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
        _glob_data(path = path, isolates= isolates, datatype = 'reads' if reads else 'contigs')
    else:
        LOGGER.critical(f"There seems to be a problem with your read path. Please provide a valid path to the reads you wish to include")
        raise SystemExit
