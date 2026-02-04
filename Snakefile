"""
SPICE Snakemake workflow

Run locally:
  SPICE_CONFIG=/path/to/config.yaml snakemake -s Snakefile --verbose -c1 --rerun-triggers mtime

Run on Slurm:
  SPICE_CONFIG=/path/to/config.yaml snakemake -s Snakefile --rerun-triggers mtime --verbose --slurm -j250 --keep-going
"""

import os
import sys
import yaml
import collections.abc
from pathlib import Path
from glob import glob

# -----------------------------
# Config loading (no spice imports)
# -----------------------------

# Snakemake sets the working directory to the Snakefile location
# Use os.getcwd() or derive from the Snakefile's actual location
# For safety, try to get the script dir, fallback to cwd
try:
    import inspect
    snakefile_dir = Path(inspect.getfile(inspect.currentframe())).resolve().parent
except Exception:
    snakefile_dir = Path(os.getcwd())

BASE_DIR = snakefile_dir
DEFAULT_CONFIG_FILE = BASE_DIR / "default_config.yaml"


def update_nested_dict(orig_dict, updated_values):
    """Update nested dictionary with overrides."""
    orig_dict = orig_dict.copy()
    for k, v in updated_values.items():
        if isinstance(v, collections.abc.Mapping):
            orig_dict[k] = update_nested_dict(orig_dict.get(k, {}), v)
        else:
            orig_dict[k] = v
    return orig_dict


def read_yaml(path):
    with open(path, "rt") as f:
        return yaml.safe_load(f.read())


snakemake_cfg = dict(config) if "config" in globals() else {}
config_path = os.environ.get("SPICE_CONFIG") or snakemake_cfg.pop("config_path", None)
selected_ids_raw = snakemake_cfg.pop("selected_ids", None)
skip_preprocessing = bool(snakemake_cfg.pop("skip_preprocessing", False))

if config_path is not None:
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file not found: {config_path}")
    user_cfg = read_yaml(config_path) or {}
elif snakemake_cfg:
    # If --configfile was used, Snakemake already populated `config`
    user_cfg = snakemake_cfg
else:
    user_cfg = {}

if not DEFAULT_CONFIG_FILE.exists():
    raise FileNotFoundError(f"Default config not found at {DEFAULT_CONFIG_FILE}")

base_cfg = read_yaml(DEFAULT_CONFIG_FILE) or {}
config_merged = update_nested_dict(base_cfg, user_cfg)
if not skip_preprocessing:
    skip_preprocessing = bool(config_merged.get("params", {}).get("skip_preprocessing", False))

# Normalize directories (similar to spice.load_config)
directories = config_merged.get("directories", {})

fallback_dirs = {
    "data_dir": BASE_DIR / "data",
    "results_dir": BASE_DIR / "results",
    "log_dir": BASE_DIR / "logs",
    "plot_dir": BASE_DIR / "plots",
    "tmp_dir": BASE_DIR / "tmp",
}

directories["base_dir"] = str(BASE_DIR)
for key, fallback in fallback_dirs.items():
    if not directories.get(key):
        directories[key] = str(fallback)

for dkey, dval in list(directories.items()):
    if isinstance(dval, str) and dval and not os.path.isabs(dval):
        directories[dkey] = str((BASE_DIR / dval).resolve())

config_merged["directories"] = directories

if "input_files" in config_merged and "knn_train" in config_merged["input_files"]:
    knn_path = config_merged["input_files"]["knn_train"]
    if knn_path and not os.path.isabs(knn_path):
        config_merged["input_files"]["knn_train"] = str((BASE_DIR / knn_path).resolve())

name = config_merged.get("name")
if not name:
    raise ValueError("Config file must specify a 'name' field.")

results_dir = os.path.join(directories["results_dir"], name)
log_dir = directories["log_dir"]

os.makedirs(results_dir, exist_ok=True)
os.makedirs(log_dir, exist_ok=True)

# Persist merged config for provenance
merged_config_path = os.path.join(results_dir, "config.yaml")
with open(merged_config_path, "wt") as f:
    yaml.safe_dump(config_merged, f)

# Prefer explicit config file for SPICE worker processes, fallback to merged
if config_path is not None:
    spice_config_path = os.path.abspath(config_path)
else:
    spice_config_path = os.path.abspath(merged_config_path)

# Optional selected IDs passed via snakemake config
selected_ids = None
if selected_ids_raw:
    if isinstance(selected_ids_raw, str):
        selected_ids = [x.strip() for x in selected_ids_raw.split(",") if x.strip()]
    elif isinstance(selected_ids_raw, (list, tuple)):
        selected_ids = [str(x) for x in selected_ids_raw]


def resolve_chrom_segments_file():
    data_dir = directories["data_dir"]
    orig = config_merged["input_files"]["copynumber"]
    processed = os.path.join(data_dir, f"{name}_processed_split.tsv")
    split_orig = orig.replace(".tsv", "_split.tsv")

    cur_file = split_orig
    if processed and os.path.exists(processed):
        cur_file = processed
    if not os.path.isabs(cur_file):
        cur_file = os.path.join(directories["base_dir"], cur_file)
    return cur_file


# -----------------------------
# Helper functions for wildcards
# -----------------------------

def get_chrom_ids_from_dir(directory):
    """Get chromosome IDs from a directory of pickle files."""
    if not os.path.exists(directory):
        return []
    return [os.path.basename(f).replace('.pickle', '') for f in glob(os.path.join(directory, '*.pickle'))]


def get_full_chrom_ids(wildcards):
    """Get all chromosome IDs from chrom_data_full directories."""
    chrom_ids = []
    for wgd in ['nowgd', 'wgd']:
        chrom_dir = checkpoints.split.get(**wildcards).output[0].replace('split.done', f'{wgd}/chrom_data_full')
        ids = get_chrom_ids_from_dir(chrom_dir)
        if selected_ids:
            ids = [x for x in ids if x in selected_ids]
        chrom_ids.extend([(wgd, cid) for cid in ids])
    return chrom_ids


def get_large_chrom_ids(wildcards):
    """Get all chromosome IDs from chrom_data_large directories."""
    chrom_ids = []
    for wgd in ['nowgd', 'wgd']:
        chrom_dir = checkpoints.split.get(**wildcards).output[0].replace('split.done', f'{wgd}/chrom_data_large')
        ids = get_chrom_ids_from_dir(chrom_dir)
        if selected_ids:
            ids = [x for x in ids if x in selected_ids]
        chrom_ids.extend([(wgd, cid) for cid in ids])
    return chrom_ids


def get_full_paths_outputs(wildcards):
    """Get all full_paths marker files that should be created."""
    outputs = []
    for wgd, cid in get_full_chrom_ids(wildcards):
        outputs.append(os.path.join(results_dir, wgd, "full_paths_started", f"{cid}.done"))
    return outputs


def get_mcmc_large_outputs(wildcards):
    """Get all mcmc_large output files that should be created."""
    outputs = []
    for wgd, cid in get_large_chrom_ids(wildcards):
        outputs.append(os.path.join(results_dir, wgd, "mcmc_solved_chroms_large", f"{cid}.pickle"))
    return outputs


# -----------------------------
# Rules (per-file with wildcards)
# -----------------------------

rule all:
    input:
        os.path.join(results_dir, "final_events.tsv"),
        os.path.join(results_dir, "summary.tsv"),


rule aggregate_full_paths:
    """Aggregate rule to ensure all full_paths jobs complete."""
    input:
        get_full_paths_outputs
    output:
        touch(os.path.join(results_dir, "full_paths_aggregate.done"))


def get_disambiguate_knn_outputs(wildcards):
    """Get knn outputs for all chromosomes with multiple solutions.
    
    This queries the checkpoint to ensure DAG is rebuilt after split,
    then scans for multiple solutions files that were created by full_paths.
    """
    # Query checkpoint to ensure we're evaluated after split completes
    checkpoints.split.get(**wildcards)
    
    outputs = []
    for wgd in ['nowgd', 'wgd']:
        multi_sol_dir = os.path.join(results_dir, wgd, 'full_paths_multiple_solutions')
        if os.path.exists(multi_sol_dir):
            for f in glob(os.path.join(multi_sol_dir, '*.pickle')):
                cid = os.path.basename(f).replace('.pickle', '')
                if selected_ids is None or cid in selected_ids:
                    outputs.append(os.path.join(results_dir, wgd, 'knn_solved_chroms', f"{cid}.pickle"))
    return outputs



rule aggregate_mcmc_large:
    """Aggregate rule to ensure all mcmc_large jobs complete."""
    input:
        get_mcmc_large_outputs
    output:
        touch(os.path.join(results_dir, "mcmc_large_aggregate.done"))


rule preprocessing:
    output:
        touch(os.path.join(results_dir, "preprocessing.done"))
    log:
        os.path.join(log_dir, "snakemake_preprocessing.log")
    resources:
        runtime=120,
        mem_mb=4000,
        job_name="spice_preprocessing"
    shell:
        r"""
        set -euo pipefail;
        if [ "{skip_preprocessing}" = "True" ]; then
            echo "Skipping preprocessing (skip_preprocessing=True)" >> {log};
            touch {output};
        else
            export SPICE_CONFIG="{spice_config_path}";
            python - <<'PY' 1>> {log} 2>> {log}
import os
import spice
spice.load_config(os.environ.get("SPICE_CONFIG"))
from spice.preprocessing.extra_preprocessing import main
from spice import config

main(
    unique_chroms=False,
    total_cn=config['params'].get('total_cn', False),
    skip_phasing=False,
    skip_centromeres=False,
)
PY
            touch {output};
        fi
        """


checkpoint split:
    input:
        rules.preprocessing.output
    output:
        touch(os.path.join(results_dir, "split.done"))
    log:
        os.path.join(log_dir, "snakemake_split.log")
    resources:
        runtime=120,
        mem_mb=4000,
        job_name="spice_split"
    shell:
        r"""
        set -euo pipefail;
        export SPICE_CONFIG="{spice_config_path}";
        python - <<'PY' 1>> {log} 2>> {log}
import os
import spice
spice.load_config(os.environ.get("SPICE_CONFIG"))
from spice.preprocessing.split_input import split_tsv_file
from spice import config

selected_ids = {selected_ids}
if isinstance(selected_ids, list) and len(selected_ids) == 0:
    selected_ids = None

split_tsv_file(
    config['name'],
    keep_old=False,
    cores=1,
    selected_ids=selected_ids,
)
PY
        touch {output}
        """


rule full_paths:
    input:
        chrom_data=os.path.join(results_dir, "{wgd}", "chrom_data_full", "{chrom_id}.pickle"),
        split_done=os.path.join(results_dir, "split.done")
    output:
        # Output will be determined dynamically (single vs multiple solutions)
        touch(os.path.join(results_dir, "{wgd}", "full_paths_started", "{chrom_id}.done"))
    log:
        os.path.join(log_dir, "full_paths_{wgd}_{chrom_id}.log")
    threads: 1
    resources:
        runtime=720,
        mem_mb=lambda wildcards: config_merged["params"].get("full_path_mem_large", 16000),
        disk_mb=lambda wildcards: config_merged["params"].get("full_path_disk", 16000),
        job_name=lambda wildcards: f"spice_full_paths_{wildcards.wgd}_{wildcards.chrom_id}"
    run:
        import os
        import sys
        import traceback
        from pathlib import Path
        sys.path.insert(0, str(BASE_DIR))
        
        try:
            os.environ["SPICE_CONFIG"] = spice_config_path
            import spice
            spice.load_config(os.environ.get("SPICE_CONFIG"))
            from spice import config
            from spice.data_loaders import resolve_data_file
            from spice.utils import timeout
            from spice.event_inference.pipeline import full_paths_from_graph_with_sv_wrapper
            from spice.logging import configure_logging, get_logger
            
            configure_logging(
                log_mode="file",
                log_dir=config['directories']['log_dir'],
                config_name=config['name'],
                level=config['params'].get('logging_level', 'INFO'),
            )
            logger = get_logger('SPICE', spice_prefix=False)
            
            chrom_segments_file = resolve_data_file()
            is_wgd = (wildcards.wgd == 'wgd')
            
            # Temporarily disable timeout to debug
            # @timeout(config['params'].get('time_limit_all_solutions', None), mode="auto")
            def run_full_paths():
                return full_paths_from_graph_with_sv_wrapper(
                    cur_id=wildcards.chrom_id,
                    is_wgd=is_wgd,
                    chrom_segments_file=chrom_segments_file,
                    sv_data_file=config['input_files'].get('sv', None),
                    chrom_file=input.chrom_data,
                    sv_matching_threshold=config['params']['sv_matching_threshold'],
                    use_cache=config['params']['use_cache'],
                    total_cn=config['params'].get('total_cn', False),
                    all_loh_solutions=config['params']['all_loh_solutions'],
                    save_output=True,
                    skip_loh_checks=True,
                )
            
            result = run_full_paths()
            logger.info(f"Successfully computed full paths for {{wildcards.chrom_id}} ({{wildcards.wgd}})")
        except Exception as e:
            print(f"ERROR in full_paths rule: {{e}}", file=sys.stderr)
            print(traceback.format_exc(), file=sys.stderr)
            # Create marker file for MCMC fallback
            fallback_dir = os.path.join(results_dir, wildcards.wgd, "full_paths_failed")
            os.makedirs(fallback_dir, exist_ok=True)
            Path(os.path.join(fallback_dir, f"{{wildcards.chrom_id}}.failed")).touch()
            raise


rule disambiguate_knn:
    input:
        full_paths=os.path.join(results_dir, "{wgd}", "full_paths_multiple_solutions", "{chrom_id}.pickle")
    output:
        os.path.join(results_dir, "{wgd}", "knn_solved_chroms", "{chrom_id}.pickle")
    log:
        os.path.join(log_dir, "disambiguate_{wgd}_{chrom_id}.log")
    threads: 1
    resources:
        runtime=720,
        mem_mb=lambda wildcards: config_merged["params"].get("knn_mem", 16000),
        job_name=lambda wildcards: f"spice_disambiguate_{wildcards.wgd}_{wildcards.chrom_id}"
    run:
        import os
        import sys
        sys.path.insert(0, str(BASE_DIR))
        
        os.environ["SPICE_CONFIG"] = spice_config_path
        import spice
        spice.load_config(os.environ.get("SPICE_CONFIG"))
        from spice import config
        from spice.data_loaders import resolve_data_file
        from spice.event_inference.pipeline import solve_with_knn_wrapper
        from spice.logging import configure_logging, get_logger
        
        configure_logging(
            log_mode="file",
            log_dir=config['directories']['log_dir'],
            config_name=config['name'],
            level=config['params'].get('logging_level', 'INFO'),
        )
        logger = get_logger('SPICE', spice_prefix=False)
        
        chrom_segments_file = resolve_data_file()
        is_wgd = (wildcards.wgd == 'wgd')
        
        full_paths_multiple_solutions_dirs = [
            os.path.join(results_dir, 'nowgd', 'full_paths_multiple_solutions'),
            os.path.join(results_dir, 'wgd', 'full_paths_multiple_solutions')
        ]
        
        result = solve_with_knn_wrapper(
            output_file=output[0],
            cur_id=wildcards.chrom_id,
            is_wgd=is_wgd,
            chrom_segments_file=chrom_segments_file,
            knn_train_data=None,
            k=config['params']['knn_k'],
            full_paths_multiple_solutions_dirs=full_paths_multiple_solutions_dirs,
            save_all_scores=None,
            perform_loh_checks=True,
            single_width_bin=True,
        )
        
        logger.info(f"Successfully disambiguated {wildcards.chrom_id} ({wildcards.wgd})")


rule mcmc_large:
    input:
        chrom_data=os.path.join(results_dir, "{wgd}", "chrom_data_large", "{chrom_id}.pickle")
    output:
        os.path.join(results_dir, "{wgd}", "mcmc_solved_chroms_large", "{chrom_id}.pickle")
    log:
        os.path.join(log_dir, "mcmc_large_{wgd}_{chrom_id}.log")
    threads: 1
    resources:
        runtime=720,
        mem_mb=8000,
        job_name=lambda wildcards: f"spice_mcmc_large_{wildcards.wgd}_{wildcards.chrom_id}"
    run:
        import os
        import sys
        sys.path.insert(0, str(BASE_DIR))
        
        os.environ["SPICE_CONFIG"] = spice_config_path
        import spice
        spice.load_config(os.environ.get("SPICE_CONFIG"))
        from spice import config
        from spice.data_loaders import resolve_data_file
        from spice.utils import timeout, FunctionTimeoutError
        from spice.event_inference.pipeline import solve_with_mcmc_wrapper
        from spice.logging import configure_logging, get_logger
        
        configure_logging(
            log_mode="file",
            log_dir=config['directories']['log_dir'],
            config_name=config['name'],
            level=config['params'].get('logging_level', 'INFO'),
        )
        logger = get_logger('SPICE', spice_prefix=False)
        
        chrom_segments_file = resolve_data_file()
        is_wgd = (wildcards.wgd == 'wgd')
        
        # Timeout handled by Snakemake's runtime resource
        # @timeout(config['params'].get('time_limit_mcmc', None), mode="auto")
        def _solve_with_mcmc_wrapper(skip_loh_check):
            return solve_with_mcmc_wrapper(
                output_file=output[0],
                chrom_file=input.chrom_data,
                is_wgd=is_wgd,
                chrom_segments_file=chrom_segments_file,
                sv_data_file=config['input_files'].get('sv', None),
                knn_train_data=None,
                k=config['params']['knn_k'],
                total_cn=config['params'].get('total_cn', False),
                save_all_scores=None,
                n_iteration_scale=config['params']['mcmc_n_iterations_scale'],
                log_progress=True,
                fail_on_empty=False,
                skip_loh_check=skip_loh_check,
            )
        
        try:
            result = _solve_with_mcmc_wrapper(skip_loh_check=False)
            logger.info(f"Successfully solved {wildcards.chrom_id} ({wildcards.wgd}) with MCMC")
        except FunctionTimeoutError:
            logger.warning(
                f'MCMC solving for {wildcards.chrom_id} timed out after {config["params"].get("time_limit_mcmc", None)} seconds. '
                'Will rerun with LOH checks skipped. Consider increasing "time_limit_mcmc" in the config file.'
            )
            result = _solve_with_mcmc_wrapper(skip_loh_check=True)
            logger.info(f"Successfully solved {wildcards.chrom_id} ({wildcards.wgd}) with MCMC (LOH checks skipped)")


rule mcmc_fallback:
    input:
        failed_marker=os.path.join(results_dir, "{wgd}", "full_paths_failed", "{chrom_id}.failed"),
        chrom_data=os.path.join(results_dir, "{wgd}", "chrom_data_full", "{chrom_id}.pickle")
    output:
        os.path.join(results_dir, "{wgd}", "mcmc_solved_chroms_full", "{chrom_id}.pickle")
    log:
        os.path.join(log_dir, "mcmc_fallback_{wgd}_{chrom_id}.log")
    threads: 1
    resources:
        runtime=720,
        mem_mb=8000,
        job_name=lambda wildcards: f"spice_mcmc_fallback_{wildcards.wgd}_{wildcards.chrom_id}"
    run:
        import os
        import sys
        sys.path.insert(0, str(BASE_DIR))
        
        os.environ["SPICE_CONFIG"] = spice_config_path
        import spice
        spice.load_config(os.environ.get("SPICE_CONFIG"))
        from spice import config
        from spice.data_loaders import resolve_data_file
        from spice.utils import timeout, FunctionTimeoutError
        from spice.event_inference.pipeline import solve_with_mcmc_wrapper
        from spice.logging import configure_logging, get_logger
        
        configure_logging(
            log_mode="file",
            log_dir=config['directories']['log_dir'],
            config_name=config['name'],
            level=config['params'].get('logging_level', 'INFO'),
        )
        logger = get_logger('SPICE', spice_prefix=False)
        logger.info(f"Running MCMC fallback for {{wildcards.chrom_id}} ({{wildcards.wgd}}) due to full_paths timeout")
        
        chrom_segments_file = resolve_data_file()
        is_wgd = (wildcards.wgd == 'wgd')
        
        # Timeout handled by Snakemake's runtime resource
        # @timeout(config['params'].get('time_limit_mcmc', None), mode="auto")
        def _solve_with_mcmc_wrapper(skip_loh_check):
            return solve_with_mcmc_wrapper(
                output_file=output[0],
                chrom_file=input.chrom_data,
                is_wgd=is_wgd,
                chrom_segments_file=chrom_segments_file,
                sv_data_file=config['input_files'].get('sv', None),
                knn_train_data=None,
                k=config['params']['knn_k'],
                total_cn=config['params'].get('total_cn', False),
                save_all_scores=None,
                n_iteration_scale=config['params']['mcmc_n_iterations_scale'],
                log_progress=True,
                fail_on_empty=False,
                skip_loh_check=skip_loh_check,
            )
        
        try:
            result = _solve_with_mcmc_wrapper(skip_loh_check=False)
            logger.info(f"Successfully solved {wildcards.chrom_id} ({wildcards.wgd}) with MCMC fallback")
        except FunctionTimeoutError:
            logger.warning(
                f'MCMC fallback for {wildcards.chrom_id} timed out after {config["params"].get("time_limit_mcmc", None)} seconds. '
                'Will rerun with LOH checks skipped.'
            )
            result = _solve_with_mcmc_wrapper(skip_loh_check=True)
            logger.info(f"Successfully solved {wildcards.chrom_id} ({wildcards.wgd}) with MCMC fallback (LOH checks skipped)")


rule combine:
    input:
        full_paths_done=os.path.join(results_dir, "full_paths_aggregate.done"),
        mcmc_large_done=os.path.join(results_dir, "mcmc_large_aggregate.done"),
        knn_outputs=get_disambiguate_knn_outputs,
    output:
        os.path.join(results_dir, "final_events.tsv"),
        os.path.join(results_dir, "summary.tsv"),
    log:
        os.path.join(log_dir, "snakemake_combine.log")
    resources:
        runtime=120,
        mem_mb=32000,
        job_name="spice_combine"
    run:
        import os
        import sys
        sys.path.insert(0, str(BASE_DIR))
        
        os.environ["SPICE_CONFIG"] = spice_config_path
        import spice
        spice.load_config(os.environ.get("SPICE_CONFIG"))
        from spice import config
        from spice.data_loaders import resolve_data_file
        from spice.event_inference.pipeline import combine_final_events
        from spice.logging import configure_logging, get_logger
        
        configure_logging(
            log_mode="file",
            log_dir=config['directories']['log_dir'],
            config_name=config['name'],
            level=config['params'].get('logging_level', 'INFO'),
        )
        logger = get_logger('SPICE', spice_prefix=False)
        
        chrom_segments_file = resolve_data_file()
        
        solved_dirs = (
            [os.path.join(results_dir, wgd, 'knn_solved_chroms') for wgd in ['nowgd', 'wgd']] +
            [os.path.join(results_dir, wgd, 'full_paths_single_solution') for wgd in ['nowgd', 'wgd']] +
            [os.path.join(results_dir, wgd, 'mcmc_solved_chroms_full') for wgd in ['nowgd', 'wgd']] +
            [os.path.join(results_dir, wgd, 'mcmc_solved_chroms_large') for wgd in ['nowgd', 'wgd']]
        )
        
        combine_final_events(
            solved_dirs=solved_dirs,
            chrom_segments_file=chrom_segments_file,
            sv_data=config['input_files'].get('sv', None),
            sv_matching_threshold=config['params']['sv_matching_threshold'],
            knn_train_data=None,
            knn_k=config['params']['knn_k'],
            output_dir=results_dir,
        )
        
        logger.info("Successfully combined all results")
