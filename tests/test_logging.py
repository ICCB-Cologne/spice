#!/usr/bin/env python
"""Test suite to verify logging setup works correctly with pytest."""

import os
import tempfile
import pytest
import logging

from spice.utils import configure_logging, get_logger, _LOGGING_CONFIG


@pytest.fixture
def temp_log_dir():
    """Create a temporary directory for log files."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield tmpdir


@pytest.fixture(autouse=True)
def reset_logging_config():
    """Reset logging configuration before each test."""
    global _LOGGING_CONFIG
    # Store original config
    original_config = _LOGGING_CONFIG.copy()
    
    # Reset to defaults
    _LOGGING_CONFIG['mode'] = 'terminal'
    _LOGGING_CONFIG['log_dir'] = None
    _LOGGING_CONFIG['config_name'] = None
    _LOGGING_CONFIG['level'] = 'INFO'
    _LOGGING_CONFIG['file_handler'] = None
    _LOGGING_CONFIG['configured'] = False
    
    yield
    
    # Restore original config
    _LOGGING_CONFIG.update(original_config)


class TestLogging:
    """Test suite for SPICE logging functionality."""
    
    def test_terminal_mode(self, temp_log_dir):
        """Test logging in terminal-only mode."""
        configure_logging(
            log_mode='terminal',
            log_dir=temp_log_dir,
            config_name='test_config',
            level='INFO'
        )
        
        logger = get_logger('TEST_TERMINAL')
        assert logger is not None
        assert logger.level == logging.INFO
        assert len(logger.handlers) == 1
        assert isinstance(logger.handlers[0], logging.StreamHandler)
        
        # Verify no log files were created
        log_files = [f for f in os.listdir(temp_log_dir) if f.endswith('.log')]
        assert len(log_files) == 0
    
    def test_file_mode(self, temp_log_dir):
        """Test logging in file-only mode."""
        configure_logging(
            log_mode='file',
            log_dir=temp_log_dir,
            config_name='test_config',
            level='INFO'
        )
        
        logger = get_logger('TEST_FILE')
        logger.info("Test message for file mode")
        
        assert logger is not None
        assert len(logger.handlers) == 1
        assert isinstance(logger.handlers[0], logging.FileHandler)
        
        # Verify log file was created
        log_files = [f for f in os.listdir(temp_log_dir) if f.endswith('.log')]
        assert len(log_files) == 1
        assert log_files[0].startswith('test_config_')
        
        # Verify log file content
        log_path = os.path.join(temp_log_dir, log_files[0])
        with open(log_path, 'r') as f:
            content = f.read()
            assert 'Test message for file mode' in content
            assert 'TEST_FILE' in content
    
    def test_both_mode(self, temp_log_dir):
        """Test logging in both terminal and file mode."""
        configure_logging(
            log_mode='both',
            log_dir=temp_log_dir,
            config_name='test_config',
            level='INFO'
        )
        
        logger = get_logger('TEST_BOTH')
        logger.info("Test message for both mode")
        
        assert logger is not None
        assert len(logger.handlers) == 2
        
        # Check we have both handler types
        handler_types = [type(h) for h in logger.handlers]
        assert logging.StreamHandler in handler_types
        assert logging.FileHandler in handler_types
        
        # Verify log file was created and contains message
        log_files = [f for f in os.listdir(temp_log_dir) if f.endswith('.log')]
        assert len(log_files) == 1
        
        log_path = os.path.join(temp_log_dir, log_files[0])
        with open(log_path, 'r') as f:
            content = f.read()
            assert 'Test message for both mode' in content
    
    def test_debug_level(self, temp_log_dir):
        """Test that DEBUG logging level is respected."""
        configure_logging(
            log_mode='terminal',
            log_dir=temp_log_dir,
            config_name='test_config',
            level='DEBUG'
        )
        
        logger = get_logger('TEST_DEBUG')
        assert logger.level == logging.DEBUG
    
    def test_multiple_loggers_share_file_handler(self, temp_log_dir):
        """Test that multiple loggers share the same file handler."""
        configure_logging(
            log_mode='file',
            log_dir=temp_log_dir,
            config_name='test_config',
            level='INFO'
        )
        
        logger1 = get_logger('LOGGER_1')
        logger2 = get_logger('LOGGER_2')
        
        logger1.info("Message from logger 1")
        logger2.info("Message from logger 2")
        
        # Only one log file should be created
        log_files = [f for f in os.listdir(temp_log_dir) if f.endswith('.log')]
        assert len(log_files) == 1
        
        # Both messages should be in the file
        log_path = os.path.join(temp_log_dir, log_files[0])
        with open(log_path, 'r') as f:
            content = f.read()
            assert 'Message from logger 1' in content
            assert 'Message from logger 2' in content
    
    def test_log_filename_format(self, temp_log_dir):
        """Test that log filename follows expected format."""
        configure_logging(
            log_mode='file',
            log_dir=temp_log_dir,
            config_name='my_test_config',
            level='INFO'
        )
        
        logger = get_logger('TEST')
        logger.info("Test message")
        
        log_files = [f for f in os.listdir(temp_log_dir) if f.endswith('.log')]
        assert len(log_files) == 1
        
        # Filename should start with config name and contain timestamp
        filename = log_files[0]
        assert filename.startswith('my_test_config_')
        assert filename.endswith('.log')
    
    def test_unconfigured_logging_defaults_to_terminal(self):
        """Test that unconfigured logging defaults to terminal mode."""
        # Don't call configure_logging
        logger = get_logger('TEST_DEFAULT', load_config=False)
        
        assert logger is not None
        assert len(logger.handlers) == 1
        assert isinstance(logger.handlers[0], logging.StreamHandler)
