import sys, pathlib, pandas, pytest, numpy,logging

from unittest.mock import patch

from bohra.SnpDetection import RunSnpDetection
# from bohra.bohra.ReRunSnpDetection import ReRunSnpDetection


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
