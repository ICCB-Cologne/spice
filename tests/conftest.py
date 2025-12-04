#!/usr/bin/env python
"""Pytest configuration and shared fixtures for SPICE tests."""

import os
import sys
import pytest

# Add the project root to the Python path
repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, repo_root)


def pytest_configure(config):
    """Configure pytest with custom markers."""
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')"
    )
    config.addinivalue_line(
        "markers", "integration: marks tests as integration tests"
    )


@pytest.fixture(scope="session")
def repo_root_dir():
    """Return the repository root directory."""
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@pytest.fixture(scope="session")
def example_data_exists(repo_root_dir):
    """Check if example data files exist."""
    example_data = os.path.join(repo_root_dir, 'data', 'example_data.tsv')
    return os.path.exists(example_data)


@pytest.fixture(scope="session")
def knn_train_data_exists(repo_root_dir):
    """Check if KNN training data exists."""
    knn_data = os.path.join(repo_root_dir, 'objects', 'knn_train_data.pickle')
    return os.path.exists(knn_data)
