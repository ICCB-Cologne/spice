#!/usr/bin/env python
"""Test suite for SPICE command-line interface."""

import os
import sys
import subprocess
import tempfile
import shutil
import pytest
import yaml


@pytest.fixture
def temp_workspace():
    """Create a temporary workspace with necessary files."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create necessary subdirectories
        data_dir = os.path.join(tmpdir, 'data')
        results_dir = os.path.join(tmpdir, 'results')
        logs_dir = os.path.join(tmpdir, 'logs')
        os.makedirs(data_dir)
        os.makedirs(results_dir)
        os.makedirs(logs_dir)
        
        # Copy example data file
        repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        example_data = os.path.join(repo_root, 'data', 'example_data.tsv')
        if os.path.exists(example_data):
            shutil.copy(example_data, os.path.join(data_dir, 'example_data.tsv'))
        
        # Create a minimal config file
        config = {
            'name': 'test_run',
            'input_files': {
                'copynumber': os.path.join(data_dir, 'example_data.tsv'),
                'knn_train': os.path.join(repo_root, 'objects', 'knn_train_data.pickle')
            },
            'directories': {
                'results_dir': results_dir,
                'log_dir': logs_dir
            },
            'params': {
                'dist_limit': 40,
                'full_path_dist_limit': 9,
                'knn_k': 250,
                'sv_matching_threshold': 10,
                'time_limit_full_paths': 60,
                'time_limit_loh_filters': 60,
                'all_loh_solutions': False,
                'use_cache': True,
                'logging_level': 'INFO',
                'mcmc_n_iterations_scale': 1.0
            }
        }
        config_path = os.path.join(tmpdir, 'test_config.yaml')
        with open(config_path, 'w') as f:
            yaml.safe_dump(config, f)
        
        yield tmpdir, config_path


class TestCLI:
    """Test suite for SPICE command-line interface."""
    
    def test_spice_command_exists(self):
        """Test that the spice command is available."""
        result = subprocess.run(
            ['spice', '--help'],
            capture_output=True,
            text=True
        )
        # Should either succeed or fail with meaningful error (not command not found)
        assert result.returncode in [0, 2]  # 0 = success, 2 = argparse error
        assert 'SPICE' in result.stdout or 'spice' in result.stdout.lower()
    
    def test_help_flag(self):
        """Test that --help flag works."""
        result = subprocess.run(
            ['spice', '--help'],
            capture_output=True,
            text=True
        )
        assert result.returncode == 0
        assert 'SPICE' in result.stdout
        assert '--config' in result.stdout
        assert '--cores' in result.stdout
        assert '--log' in result.stdout
    
    def test_config_required(self):
        """Test that config argument is required."""
        result = subprocess.run(
            ['spice'],
            capture_output=True,
            text=True
        )
        # Should fail without config
        assert result.returncode != 0
        assert '--config' in result.stderr or 'required' in result.stderr.lower()
    
    def test_config_file_not_found(self):
        """Test that missing config file produces appropriate error."""
        result = subprocess.run(
            ['spice', '--config', '/nonexistent/config.yaml'],
            capture_output=True,
            text=True
        )
        assert result.returncode != 0

    def test_invalid_step_rejected(self, temp_workspace):
        """Test that invalid step names are rejected."""
        tmpdir, config_path = temp_workspace
        
        result = subprocess.run(
            ['spice', 'invalid_step', '--config', config_path],
            capture_output=True,
            text=True
        )
        assert result.returncode != 0
        assert 'invalid' in result.stderr.lower() or 'Invalid' in result.stderr
    
    def test_valid_steps_accepted(self, temp_workspace):
        """Test that valid step names are accepted."""
        tmpdir, config_path = temp_workspace
        
        valid_steps = ['split', 'all_solutions', 'disambiguate', 'large_chroms', 'combine']
        for step in valid_steps:
            # Just test that the command parses correctly (will fail on execution)
            # Use --clean to avoid actually running
            result = subprocess.run(
                ['spice', step, '--config', config_path, '--clean'],
                capture_output=True,
                text=True,
                timeout=10
            )
            # Should not fail on argument parsing
            assert 'Invalid step' not in result.stderr
    
    def test_valid_steps_with_plus_accepted(self, temp_workspace):
        """Test that valid step names are accepted."""
        tmpdir, config_path = temp_workspace
        
        valid_steps = ['split', 'all_solutions', 'disambiguate', 'large_chroms', 'combine']
        for step in valid_steps:
            # Just test that the command parses correctly (will fail on execution)
            # Use --clean to avoid actually running
            result = subprocess.run(
                ['spice', step + '+', '--config', config_path, '--clean'],
                capture_output=True,
                text=True,
                timeout=10
            )
            # Should not fail on argument parsing
            assert 'Invalid step' not in result.stderr
    
    def test_all_step_works(self, temp_workspace):
        """Test that 'all' step is accepted."""
        tmpdir, config_path = temp_workspace
        
        result = subprocess.run(
            ['spice', 'all', '--config', config_path, '--clean'],
            capture_output=True,
            text=True,
            timeout=10
        )
        # Should not fail on argument parsing
        assert 'Invalid step' not in result.stderr
    
    def test_cores_flag(self, temp_workspace):
        """Test that --cores flag is accepted."""
        tmpdir, config_path = temp_workspace
        
        result = subprocess.run(
            ['spice', '--config', config_path, '--cores', '4', '--clean'],
            capture_output=True,
            text=True,
            timeout=10
        )
        # Should not fail on argument parsing
        assert 'cores' not in result.stderr.lower() or result.returncode == 0
    
    def test_log_flag_terminal(self, temp_workspace):
        """Test that --log terminal flag works."""
        tmpdir, config_path = temp_workspace
        
        result = subprocess.run(
            ['spice', '--config', config_path, '--log', 'terminal', '--clean'],
            capture_output=True,
            text=True,
            timeout=10
        )
        # Should accept the flag
        assert result.returncode == 0 or 'log' not in result.stderr.lower()
    
    def test_log_flag_file(self, temp_workspace):
        """Test that --log file flag works."""
        tmpdir, config_path = temp_workspace
        
        result = subprocess.run(
            ['spice', '--config', config_path, '--log', 'file', '--clean'],
            capture_output=True,
            text=True,
            timeout=10
        )
        # Should accept the flag
        assert result.returncode == 0 or 'log' not in result.stderr.lower()
    
    def test_log_flag_both(self, temp_workspace):
        """Test that --log both flag works."""
        tmpdir, config_path = temp_workspace
        
        result = subprocess.run(
            ['spice', '--config', config_path, '--log', 'both', '--clean'],
            capture_output=True,
            text=True,
            timeout=10
        )
        # Should accept the flag
        assert result.returncode == 0 or 'log' not in result.stderr.lower()
    
    def test_debug_flag(self, temp_workspace):
        """Test that --debug flag is accepted."""
        tmpdir, config_path = temp_workspace
        
        result = subprocess.run(
            ['spice', '--config', config_path, '--debug', '--clean'],
            capture_output=True,
            text=True,
            timeout=10
        )
        # Should not fail on argument parsing
        assert result.returncode == 0 or 'debug' not in result.stderr.lower()
    
    def test_keep_old_flag(self, temp_workspace):
        """Test that --keep-old flag is accepted."""
        tmpdir, config_path = temp_workspace
        
        result = subprocess.run(
            ['spice', '--config', config_path, '--keep-old', '--clean'],
            capture_output=True,
            text=True,
            timeout=10
        )
        # Should not fail on argument parsing
        assert result.returncode == 0 or 'keep-old' not in result.stderr.lower()
    
    def test_clean_flag(self, temp_workspace):
        """Test that --clean flag works and cleans directories."""
        tmpdir, config_path = temp_workspace
        results_dir = os.path.join(tmpdir, 'results', 'test_run')
        
        # Create some dummy directories
        os.makedirs(os.path.join(results_dir, 'wgd'), exist_ok=True)
        os.makedirs(os.path.join(results_dir, 'nowgd'), exist_ok=True)
        dummy_file = os.path.join(results_dir, 'wgd', 'test.txt')
        with open(dummy_file, 'w') as f:
            f.write('test')
        
        result = subprocess.run(
            ['spice', '--config', config_path, '--clean'],
            capture_output=True,
            text=True
        )
        
        # Should execute successfully
        assert result.returncode == 0
        assert 'Cleaning intermediate files' in result.stdout or 'Cleaning intermediate files' in result.stderr
        # Directories should be cleaned
        assert not os.path.exists(dummy_file)
    
    def test_ids_flag(self, temp_workspace):
        """Test that --ids flag is accepted."""
        tmpdir, config_path = temp_workspace
        
        result = subprocess.run(
            ['spice', '--config', config_path, '--ids', 'sample1,sample2', '--clean'],
            capture_output=True,
            text=True,
            timeout=10
        )
        # Should not fail on argument parsing
        assert result.returncode == 0 or 'ids' not in result.stderr.lower()
    
    def test_multiple_steps(self, temp_workspace):
        """Test that multiple steps can be specified."""
        tmpdir, config_path = temp_workspace
        
        result = subprocess.run(
            ['spice', 'split', 'disambiguate', '--config', config_path, '--clean'],
            capture_output=True,
            text=True,
            timeout=10
        )
        # Should not fail on argument parsing
        assert 'Invalid step' not in result.stderr
