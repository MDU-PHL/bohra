import subprocess
import pytest
from unittest.mock import patch, MagicMock
from bohra.launcher.CheckDeps import check_dependencies, _check_databases


def test_check_dependencies_stderr_not_captured():
    """
    Test that check_dependencies does not capture stderr,
    allowing stderr output from bash script to pass through.
    """
    with patch('subprocess.Popen') as mock_popen:
        # Setup mock process
        mock_process = MagicMock()
        mock_process.poll.side_effect = [None, 0]  # Running, then complete
        mock_process.stdout.readline.side_effect = ["test output\n", ""]
        mock_process.returncode = 0
        mock_popen.return_value = mock_process
        
        # Call the function
        result = check_dependencies()
        
        # Verify Popen was called without stderr=subprocess.PIPE
        call_args = mock_popen.call_args
        assert 'stderr' not in call_args[1] or call_args[1].get('stderr') is None
        assert call_args[1]['stdout'] == subprocess.PIPE
        assert result == 0


def test_check_databases_stderr_not_captured():
    """
    Test that _check_databases does not capture stderr,
    allowing stderr output from bash script to pass through.
    """
    with patch('subprocess.Popen') as mock_popen:
        # Setup mock process
        mock_process = MagicMock()
        mock_process.poll.side_effect = [None, 0]  # Running, then complete
        mock_process.stdout.readline.side_effect = ["checking dbs\n", ""]
        mock_process.returncode = 0
        mock_popen.return_value = mock_process
        
        # Call the function
        result = _check_databases(db_install=False)
        
        # Verify Popen was called without stderr=subprocess.PIPE
        call_args = mock_popen.call_args
        assert 'stderr' not in call_args[1] or call_args[1].get('stderr') is None
        assert call_args[1]['stdout'] == subprocess.PIPE
        assert result == 0


def test_check_dependencies_error_handling():
    """
    Test that check_dependencies properly handles errors
    without trying to read from uncaptured stderr.
    """
    with patch('subprocess.Popen') as mock_popen:
        # Setup mock process that fails
        mock_process = MagicMock()
        mock_process.poll.side_effect = [None, 1]  # Running, then error
        mock_process.stdout.readline.side_effect = ["error output\n", ""]
        mock_process.returncode = 1
        mock_popen.return_value = mock_process
        
        # Call the function
        result = check_dependencies()
        
        # Should return 1 on error
        assert result == 1


def test_check_databases_error_handling():
    """
    Test that _check_databases properly handles errors
    without trying to read from uncaptured stderr.
    """
    with patch('subprocess.Popen') as mock_popen:
        # Setup mock process that fails
        mock_process = MagicMock()
        mock_process.poll.side_effect = [None, 1]  # Running, then error
        mock_process.stdout.readline.side_effect = ["db error\n", ""]
        mock_process.returncode = 1
        mock_popen.return_value = mock_process
        
        # Call the function
        result = _check_databases(db_install=True)
        
        # Should return 1 on error
        assert result == 1
