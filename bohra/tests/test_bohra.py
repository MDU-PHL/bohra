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
              
#                 detect_obj = RunSnpDetection()
#                 detect_obj.logger = logging.getLogger(__name__)
#                 tab = pandas.DataFrame({'A':[1,2,3]})
#                 assert detect_obj.min_four_samples(tab)

# def test_missing():
#         '''
#         a full dataframe returns true
#         '''
#         with patch.object(RunSnpDetection, "__init__", lambda x: None):
#                 detect_obj = RunSnpDetection()
#                 detect_obj.logger = logging.getLogger(__name__)
#                 tab = pandas.DataFrame({'A':[1,2,3,4], 'B':[5,6,7,8], 'C':[9,10,11,12]})
#                 assert detect_obj.all_data_filled(tab)

# def test_with_missing():
#         '''
#         return False if missing data present
#         '''
#         with patch.object(RunSnpDetection, "__init__", lambda x: None):
#                 detect_obj = RunSnpDetection()
#                 detect_obj.logger = logging.getLogger(__name__)
#                 tab = pandas.DataFrame({'A':[1,2,4, numpy.nan], 'B':[5,6,7,8], 'C':[9,10,11,12]})
#                 assert detect_obj.all_data_filled(tab) == False


# def test_structure():
#         with patch.object(RunSnpDetection, "__init__", lambda x: None):
#                 detect_obj = RunSnpDetection()
#                 detect_obj.logger = logging.getLogger(__name__)
#                 tab = pandas.DataFrame({'A':[1,2,3,4], 'B':[5,6,7,8], 'C':[9,10,11,12]})
#                 assert detect_obj.check_input_structure(tab) == True
=======

# def test_path_exists():
#         '''
#         test that path_exists returns True
#         '''
#         with patch.object(RunSnpDetection, "__init__", lambda x: None):
#                 p = pathlib.Path('bohra','tests', 'test_file.txt')
#                 detect_obj = RunSnpDetection()
#                 detect_obj.logger = logging.getLogger(__name__)
#                 assert detect_obj.path_exists(p)

