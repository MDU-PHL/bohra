import sys, pathlib, pandas, pytest, numpy,logging, subprocess,os, json

from unittest.mock import patch

from bohra.SnpDetection import RunSnpDetection
# from bohra.bohra.ReRunSnpDetection import ReRunSnpDetection

script_path = pathlib.Path(__file__).parent.parent / 'tests'
conda_path = pathlib.Path(f"{os.environ['CONDA_PREFIX']}").parent
outdir = pathlib.Path.cwd()

with open(f"{script_path / 'variables.json'}", 'r') as j:
        cfg = json.load(j)


def get_main(tool):

        return cfg[f"{tool}"]["main"]

def get_input(tool):
        return cfg[f"{tool}"]["input"]

def get_truth(tool):
        return cfg[f"{tool}"]["truth"]

def get_paths(tool):
        nf = script_path / get_main(tool)
        cp= script_path / get_input(tool)
        tr = script_path / get_truth(tool)
        return nf,cp,tr

def clean():

        subprocess.run(f"rm -rf .nextflow* bohra.log work for_testing", shell = True)
        


@pytest.mark.nf
def test_nf_emmtyper():
        nf,cp,tr = get_paths("emmtyper")
        
        cmd = f"nextflow run {nf} --contig_path {cp} --conda_path {conda_path} --outdir {outdir} --publish_dir_mode copy --enable_conda true -with-conda --truth {tr}"
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        
        assert proc.returncode == 0
        # clean()

@pytest.mark.nf
def test_nf_ectyper():
        nf,cp,tr = get_paths("ectyper")
        
        cmd = f"nextflow run {nf} --contig_path {cp} --conda_path {conda_path} --outdir {outdir} --publish_dir_mode copy --enable_conda true -with-conda --truth {tr}"
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        
        assert proc.returncode == 0
        clean()

@pytest.mark.nf
def test_nf_mlst():
        nf,cp,tr = get_paths("mlst")
        
        cmd = f"nextflow run {nf} --contig_path {cp} --conda_path {conda_path} --outdir {outdir} --publish_dir_mode copy --enable_conda true -with-conda --blast_db no_db --data_dir no_db --mlst_exclude 'ecoli,abaumannii' --truth {tr}"
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        
        assert proc.returncode == 0
        clean()

def test_name_string():
        '''
        assert true when the input is a string > len 0
        '''
        with patch.object(RunSnpDetection, "__init__", lambda x: None):
                detect_obj = RunSnpDetection()
                detect_obj.logger = logging.getLogger(__name__)
                assert detect_obj._name_exists('somestring')

def test_name_non_string():
        '''
        if a non string input is used return false
        '''
        with patch.object(RunSnpDetection, "__init__", lambda x: None):
                detect_obj = RunSnpDetection()
                detect_obj.logger = logging.getLogger(__name__)
                with pytest.raises(SystemExit):
                        detect_obj._name_exists(9)

def test_name_empty_string():
        '''
        confirm that False is returned if a 0 length strin is input
        '''
        with patch.object(RunSnpDetection, "__init__", lambda x: None):
                detect_obj = RunSnpDetection()
                detect_obj.logger = logging.getLogger(__name__)
                with pytest.raises(SystemExit):
                        detect_obj._name_exists('')
              
def test_path_exists():

        '''
        confirm true is returned when path is found
        '''
        with patch.object(RunSnpDetection, "__init__", lambda x: None):
                detect_obj = RunSnpDetection()
                detect_obj.logger = logging.getLogger(__name__)
                some_path = pathlib.Path(__file__).parent
                assert detect_obj._path_exists(some_path)

            
def test_not_path_exists():

        '''
        confirm true is returned when path is found
        '''
        with patch.object(RunSnpDetection, "__init__", lambda x: None):
                detect_obj = RunSnpDetection()
                detect_obj.logger = logging.getLogger(__name__)
                some_path = pathlib.Path('some_file.txt')
                assert not detect_obj._path_exists(some_path)
