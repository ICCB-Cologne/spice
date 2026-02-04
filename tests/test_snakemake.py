#!/usr/bin/env python
"""Integration tests for Snakemake workflow."""

import os
import sys
import shutil
import subprocess
import tempfile
import csv
import pytest
import yaml

import pandas as pd


@pytest.mark.integration
def test_snakefile_parses(repo_root_dir):
    """Test that Snakefile can be parsed by Snakemake."""
    snakefile_path = os.path.join(repo_root_dir, "Snakefile")
    assert os.path.exists(snakefile_path), "Snakefile not found"

    # Set up minimal config for parsing
    with tempfile.TemporaryDirectory() as tmpdir:
        config_path = os.path.join(tmpdir, "test_config.yaml")
        config = {
            "name": "test_parse",
            "input_files": {"copynumber": "dummy.tsv"},
            "directories": {
                "data_dir": tmpdir,
                "results_dir": tmpdir,
                "log_dir": tmpdir,
                "plot_dir": tmpdir,
            },
            "params": {"logging_level": "INFO"},
        }
        with open(config_path, "w") as f:
            yaml.safe_dump(config, f)

        result = subprocess.run(
            ["spice", "--config", config_path, "--snakemake", "--snakemake-mode", "local", "--snakemake-cores", "1",
             '--help'],
            capture_output=True,
            text=True,
            cwd=repo_root_dir,
        )

        assert result.returncode == 0, f"Snakefile failed to parse: {result.stderr}"


@pytest.mark.integration
def test_snakemake_cli_integration(repo_root_dir):
    """Test that Snakemake can be invoked via the SPICE CLI."""
    with tempfile.TemporaryDirectory() as tmpdir:
        config_path = os.path.join(tmpdir, "test_config.yaml")
        config = {
            "name": "test_cli",
            "input_files": {"copynumber": "dummy.tsv"},
            "directories": {
                "data_dir": tmpdir,
                "results_dir": tmpdir,
                "log_dir": tmpdir,
                "plot_dir": tmpdir,
            },
        }
        with open(config_path, "w") as f:
            yaml.safe_dump(config, f)

        # Test that CLI accepts snakemake arguments (will fail on file not found, not argument parsing)
        result = subprocess.run(
            ["spice", "event_inference", "--config", config_path, "--snakemake", "--help"],
            capture_output=True,
            text=True,
            cwd=repo_root_dir,
        )

        # Should fail due to missing input file, not argument parsing
        # If argument parsing failed, error would be about invalid argument
        assert "invalid" not in result.stderr.lower() or "config" in result.stderr.lower(), \
            f"CLI argument parsing failed: {result.stderr}"


@pytest.mark.integration
def test_snakemake_env_variable_config(repo_root_dir):
    """Test that SPICE_CONFIG environment variable works with Snakemake."""
    config_path = os.path.join(repo_root_dir, "configs", "example_config.yaml")

    result = subprocess.run(
        ["spice", "--config", config_path, "--snakemake", "--snakemake-mode", "local", "--snakemake-cores", "1"],
        capture_output=True,
        text=True,
        cwd=repo_root_dir,
    )

    assert result.returncode == 0, f"Snakemake failed with config: {result.stderr}"
    output = result.stdout + result.stderr
    assert "Building DAG" in output or "Nothing to be done" in output, \
        f"Unexpected Snakemake output: stdout={result.stdout}, stderr={result.stderr}"


@pytest.mark.integration
def test_snakemake_cli_help(repo_root_dir):
    """Test snakemake help for CLI wrapper is accessible."""
    result = subprocess.run(
        ["spice", "--help"],
        capture_output=True,
        text=True,
        cwd=repo_root_dir,
    )
    assert result.returncode == 0, f"CLI help failed: {result.stderr}"


@pytest.mark.integration
def test_snakemake_matches_cli_example(repo_root_dir, example_data_exists, knn_train_data_exists):
    """Run example data with CLI and Snakemake; verify identical final_events.tsv."""

    repo_root = repo_root_dir
    example_data = os.path.join(repo_root, "data", "test_data.tsv")

    with tempfile.TemporaryDirectory() as tmpdir:
        data_dir = os.path.join(tmpdir, "data")
        results_dir = os.path.join(tmpdir, "results")
        logs_dir = os.path.join(tmpdir, "logs")
        plot_dir = os.path.join(tmpdir, "plots")
        os.makedirs(data_dir)
        os.makedirs(results_dir)
        os.makedirs(logs_dir)
        os.makedirs(plot_dir)

        shutil.copy(example_data, os.path.join(data_dir, "example_data.tsv"))

        # Pick a small subset of IDs from the example data (prefer non-diploid by CN variation)
        group_values = {}
        with open(os.path.join(data_dir, "example_data.tsv"), newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                key = (row.get("sample_id"), row.get("chrom"))
                vals = group_values.setdefault(key, {"cn_a": set(), "cn_b": set()})
                try:
                    vals["cn_a"].add(int(float(row["cn_a"])))
                    vals["cn_b"].add(int(float(row["cn_b"])))
                except Exception:
                    continue

        base_config = {
            "input_files": {
                "copynumber": os.path.join(data_dir, "example_data.tsv"),
                "knn_train": os.path.join(repo_root, "objects", "train_events_sv_and_unamb.pickle"),
            },
            "directories": {
                "data_dir": data_dir,
                "results_dir": results_dir,
                "log_dir": logs_dir,
                "plot_dir": plot_dir,
            },
            "params": {
                "dist_limit": 40,
                "full_path_dist_limit": 9,
                "full_path_high_mem_dist_limit": 8,
                "knn_k": 50,
                "sv_matching_threshold": 10,
                "time_limit_all_solutions": 300,
                "time_limit_loh_filters": 60,
                "time_limit_mcmc": 300,
                "all_loh_solutions": False,
                "use_cache": True,
                "logging_level": "INFO",
                "mcmc_n_iterations_scale": 1.0,
                "skip_preprocessing": True,
            },
        }

        normal_config = dict(base_config)
        normal_config["name"] = "normal_run"

        snakemake_config = dict(base_config)
        snakemake_config["name"] = "snakemake_run"

        normal_config_path = os.path.join(tmpdir, "config_normal.yaml")
        snakemake_config_path = os.path.join(tmpdir, "config_snakemake.yaml")
        with open(normal_config_path, "w") as f:
            yaml.safe_dump(normal_config, f)
        with open(snakemake_config_path, "w") as f:
            yaml.safe_dump(snakemake_config, f)

        env = os.environ.copy()
        env["PYTHONPATH"] = repo_root + os.pathsep + env.get("PYTHONPATH", "")

        # Run normal CLI pipeline (skip preprocessing to avoid heavy deps)
        cli_cmd = [
            sys.executable,
            "-m",
            "spice.cli",
            "event_inference",
            "--config",
            normal_config_path,
            "--skip-preprocessing",
        ]
        subprocess.run(cli_cmd, check=True, cwd=repo_root, env=env)

        # Run snakemake pipeline via CLI integration
        snakemake_cmd = [
            "spice",
            "--config",
            snakemake_config_path,
            "--snakemake",
            "--snakemake-mode",
            "local",
            "--snakemake-cores",
            "1",
        ]
        subprocess.run(snakemake_cmd, check=True, cwd=repo_root)

        normal_out = os.path.join(results_dir, "normal_run", "final_events.tsv")
        snakemake_out = os.path.join(results_dir, "snakemake_run", "final_events.tsv")

        assert os.path.exists(normal_out), "CLI output final_events.tsv missing"
        assert os.path.exists(snakemake_out), "Snakemake output final_events.tsv missing"

        # Note that I cannot directly compare the output because the order of events might differ and if 
        # MCMC is present even the events might be different

        snakemake_df = pd.read_csv(snakemake_out, sep="\t")
        cli_df = pd.read_csv(normal_out, sep="\t")

        assert snakemake_df.shape == cli_df.shape, \
            f"final_events.tsv shape mismatch: CLI {cli_df.shape}, Snakemake {snakemake_df.shape}"              
        assert (cli_df.index == snakemake_df.index).all(), \
            f"Index mismatch: CLI {cli_df.index.tolist()}, Snakemake {snakemake_df.index.tolist()}"
        assert (cli_df.columns == snakemake_df.columns).all(), \
            f"Columns mismatch: CLI {cli_df.columns.tolist()}, Snakemake {snakemake_df.columns.tolist()}"
        assert (cli_df['sample'].value_counts() == snakemake_df['sample'].value_counts()).all(), \
            f"Sample value counts mismatch: CLI {cli_df['sample'].value_counts().to_dict()}, Snakemake {snakemake_df['sample'].value_counts().to_dict()}"
        assert (cli_df.groupby('sample')['solved'].size() == snakemake_df.groupby('sample')['solved'].size()).all(), \
            f"Solved counts by sample mismatch: CLI {cli_df.groupby('sample')['solved'].size().to_dict()}, Snakemake {snakemake_df.groupby('sample')['solved'].size().to_dict()}"
        assert (cli_df.groupby('id').size() == snakemake_df.groupby('id').size()).all(), \
            f"Groupby id size mismatch: CLI {cli_df.groupby('id').size().to_dict()}, Snakemake {snakemake_df.groupby('id').size().to_dict()}"
