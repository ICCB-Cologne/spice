#!/usr/bin/env python
"""Command-line interface for SPICE."""

import os
import sys
import argparse
import subprocess
import re

# Import base package only; defer submodule imports until after config is loaded
import spice


def get_version():
    """Extract version from setup.py."""
    setup_path = os.path.join(os.path.dirname(__file__), '..', 'setup.py')
    with open(setup_path, 'r') as f:
        content = f.read()
    match = re.search(r"version\s*=\s*['\"]([^'\"]+)['\"]", content)
    if match:
        return match.group(1)
    return 'unknown'


def main_event_inference(args):
    """Run event inference pipeline."""
    # Handle 'all' or empty arguments, and expand trailing + syntax (e.g., split+)
    step_order = ['preprocessing', 'split', 'all_solutions', 'disambiguate', 'large_chroms', 'combine']
    valid_steps = step_order + ['all']
    steps_input = getattr(args, 'event_steps', None)
    
    if steps_input is None or steps_input == ['all']:
        which = step_order.copy()
    elif any('+' in x for x in steps_input):
        assert len(steps_input) == 1, 'Only one step with + is allowed'
        assert steps_input[0].endswith('+'), 'Only trailing + syntax is supported'
        which = step_order[step_order.index(steps_input[0][:-1]):]
    else:
        which = steps_input

    invalid_steps = [step for step in which if step not in valid_steps]
    if invalid_steps:
        raise ValueError(f"Invalid step(s): {', '.join(invalid_steps)}. Valid steps are: preprocessing, split, all_solutions, disambiguate, large_chroms, combine")

    # Handle unlock early to avoid expensive imports
    if args.unlock:
        spice.set_config(args.config_path)
        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
        snakefile = os.path.join(repo_root, 'Snakefile_event_inference')
        if not os.path.exists(snakefile):
            raise FileNotFoundError(f"Snakefile_event_inference not found at {snakefile}")

        cmd = [
            'snakemake',
            '-s', snakefile,
            '--configfile', args.config_path,
            '--unlock'
        ]

        env = os.environ.copy()
        env['SPICE_CONFIG'] = os.path.abspath(args.config_path)
        
        print(f"Unlocking Snakemake working directory with config: {args.config_path}")
        result = subprocess.run(cmd, env=env)
        
        if result.returncode == 0:
            print("Successfully unlocked Snakemake working directory.")
        else:
            print(f"Unlock failed with return code {result.returncode}")
            sys.exit(result.returncode)
        
        return

    # Handle snakemake mode early to avoid expensive imports
    if args.snakemake:
        spice.set_config(args.config_path)
        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
        snakefile = os.path.join(repo_root, 'Snakefile_event_inference')
        if not os.path.exists(snakefile):
            raise FileNotFoundError(f"Snakefile_event_inference not found at {snakefile}")

        if args.skip_preprocessing:
            spice.load_config(args.config_path)
            from spice import config
            import shutil

            name = config.get('name')
            if not name:
                raise ValueError("Config file must specify a 'name' field.")
            data_dir = config['directories']['data_dir']
            src = config['input_files']['copynumber']
            dst = os.path.join(data_dir, f"{name}_processed.tsv")
            os.makedirs(data_dir, exist_ok=True)

            if not os.path.isabs(src):
                src = os.path.join(config['directories']['base_dir'], src)

            if not os.path.exists(dst):
                print(f"--skip-preprocessing set. Copying {src} -> {dst}")
                shutil.copyfile(src, dst)
            else:
                print(f"--skip-preprocessing set. Using existing {dst}")

        cmd = [
            'snakemake',
            '-s', snakefile,
            '--rerun-triggers', 'mtime',
            '--verbose',
            '--configfile', args.config_path,
            '--config', f'config_path={args.config_path}',
            '--keep-going'
        ]
        
        # Pass skip_preprocessing to snakemake if set
        if args.skip_preprocessing:
            cmd.extend(['--config', 'skip_preprocessing=True'])
        
        # add execution mode and number of jobs/cores
        if args.snakemake_mode == 'slurm':
            cmd.extend(['--slurm', '-j', str(args.snakemake_jobs)])
        else:
            # Local mode: explicitly disable profiles and cluster submission
            cmd.extend(['--profile', '', '--cores', str(args.snakemake_cores)])

        env = os.environ.copy()
        env['SPICE_CONFIG'] = os.path.abspath(args.config_path)
        subprocess.run(cmd, check=True, env=env)
        return

    # Load configuration before importing submodules that may read it
    spice.load_config(args.config_path)
    from spice import config
    
    # Handle --clean early to avoid expensive imports
    if args.clean:
        import shutil
        name = config.get('name', None)
        if not name:
            raise ValueError("Config file must specify a 'name' field.")
        results_dir = os.path.join(config['directories']['results_dir'], name)
        print(f'Cleaning intermediate files at {results_dir}')
        for wgd in ['nowgd', 'wgd']:
            shutil.rmtree(os.path.join(results_dir, wgd), ignore_errors=True)
        print('Done cleaning.')
        return
    
    # Now do the expensive imports
    from spice.data_loaders import load_final_events, resolve_data_file
    from spice.utils import timeout, FunctionTimeoutError
    from spice.logging import configure_logging, get_logger
    from spice.cli_functions import save_fail_reports, step_aware_cleanup, _run_batch
    from spice.preprocessing.split_input import split_tsv_file
    from spice.event_inference.pipeline import (
        full_paths_from_graph_with_sv_wrapper, solve_with_knn_wrapper, solve_with_mcmc_wrapper,
        combine_final_events)

    if 'name' not in config or not config['name']:
        raise ValueError("Config file must specify a 'name' field.")
    if 'input_files' not in config or 'copynumber' not in config['input_files']:
        raise ValueError("Config file must specify 'input_files.copynumber'.")
    
    # Create logger AFTER imports to avoid it being disabled by medicc's logging.config.dictConfig
    log_level = 'DEBUG' if args.debug else config['params'].get('logging_level', 'INFO')

    configure_logging(
        log_mode=args.log,
        log_dir=config['directories']['log_dir'],
        config_name=config['name'],
        level=log_level,
    )
    logger = get_logger('SPICE', spice_prefix=False)

    if args.debug and args.cores and args.cores > 1:
        logger.warning("Debug mode with multiple cores may lead to interleaved log messages.")

    name = config['name']
    if ' ' in name:
        logger.error("Project name must not contain spaces.")
        return
    directories = config['directories']
    results_dir = os.path.join(directories['results_dir'], name)
    log_dir = os.path.join(directories['log_dir'])
    plots_base_dir = os.path.join(directories['plot_dir'], name)
    for cur_dir in [results_dir, log_dir, plots_base_dir]:
        if not os.path.exists(cur_dir):
            logger.info(f"Creating directory {cur_dir}")
            os.makedirs(cur_dir)

    logger.info('Running SPICE: Selection Patterns In somatic Copy-number Events')
    logger.info(f'Running event inference for project name {name} with config file {args.config_path}')

    logger.info(f'Results will be stored in {results_dir}')
    logger.info(f'Running the following steps: {", ".join(which)}')

    selected_ids = args.ids.split(',') if args.ids is not None else None
    if selected_ids is not None:
        logger.info(f'Selecting only IDs: {selected_ids}')

    # Check number of samples and warn if many
    try:
        import pandas as pd
        copynumber_file = config['input_files']['copynumber']
        df = pd.read_csv(copynumber_file, sep='\t', usecols=['sample_id'])
        n_samples = df['sample_id'].nunique()
        logger.info(f'Input copy-number file contains {n_samples} unique samples.')
        if n_samples > 50:
            logger.warning("=" * 80)
            logger.warning("!!! WARNING !!!")
            logger.warning("Large number of input samples detected (N={n_samples}).")
            logger.warning("")
            logger.warning("SPICE can be very slow when processing many samples in serial mode.")
            logger.warning("For large datasets, we strongly recommend using the Snakemake workflow")
            logger.warning("for parallel execution on a cluster.")
            logger.warning("")
            logger.warning("See README section 'Using with Snakemake' for more information:")
            logger.warning("=" * 80)
    except Exception as e:
        logger.debug(f"Could not count samples from input file: {e}")

    total_cn = config['params'].get('total_cn', False)
    if total_cn:
        raise NotImplementedError("total_cn=True is not yet supported in this version of SPICE")

    # Clean old files
    if not args.keep_old and args.ids is None:
        if config['params'].get('skip_existing', False):
            raise ValueError("If skip_existing=True in config, have to use --keep-old to avoid deleting existing files.")
        logger.info('Cleaning old intermediate files')
        step_aware_cleanup(results_dir, which)

    # Run preprocessing first unless skipped
    if 'preprocessing' in which and not args.skip_preprocessing:
        from spice.preprocessing.extra_preprocessing import main as extra_preprocessing_main
        logger.info('Starting extra preprocessing step (pre-split)')
        extra_preprocessing_main(
            unique_chroms=bool(args.pre_unique_chroms),
            total_cn=config['params'].get('total_cn', False),
            skip_phasing=bool(args.pre_skip_phasing),
            skip_centromeres=bool(args.pre_skip_centromeres),
        )
    elif 'preprocessing' in which and args.skip_preprocessing:
        logger.info('Skipping preprocessing due to --skip-preprocessing')

    chrom_segments_file = resolve_data_file()

    if 'split' in which:
        logger.info('Starting splitting of the input TSV')
        split_tsv_file(name, keep_old=args.keep_old, cores=args.cores, selected_ids=selected_ids)

    # Collect per-ID failures to report at end
    failed_reports = []

    if 'all_solutions' in which:
        logger.info('Starting inference of all solutions')
        skip_existing = config['params'].get('skip_existing', False)
        for wgd_status in ['nowgd', 'wgd']:
            is_wgd = (wgd_status == 'wgd')
            cur_ids = [x.replace('.pickle', '')
                    for x in os.listdir(os.path.join(str(results_dir), wgd_status, 'chrom_data_full'))]
            if selected_ids is not None:
                cur_ids = [x for x in cur_ids if x in selected_ids]

            @timeout(config['params']['time_limit_all_solutions'], mode="auto")
            def run_full_paths(cur_id):
                output_file = os.path.join(results_dir, wgd_status, 'full_paths_multiple_solutions', f'{cur_id}.pickle')
                if skip_existing and os.path.exists(output_file):
                    logger.info(f"Skipping all_solutions for {cur_id} ({wgd_status}) as {output_file} exists.")
                    return {'status': 'skipped', 'cur_id': cur_id, 'step': 'all_solutions'}
                return full_paths_from_graph_with_sv_wrapper(
                    cur_id=cur_id,
                    is_wgd=(wgd_status == 'wgd'),
                    chrom_segments_file=chrom_segments_file,
                    sv_data_file=config['params'].get('sv_data_file', None),
                    chrom_file=os.path.join(results_dir, wgd_status, 'chrom_data_full', f'{cur_id}.pickle'),
                    sv_matching_threshold=config['params']['sv_matching_threshold'],
                    use_cache=config['params']['use_cache'],
                    total_cn=config['params'].get('total_cn', False),
                    all_loh_solutions=config['params']['all_loh_solutions'],
                    save_output=True,
                    skip_loh_checks=True,
                )
            results = _run_batch(cur_ids, args.cores, f'All solutions ({wgd_status})', run_full_paths, logger)
            cur_failed_reports = [r for r in results if isinstance(r, dict) and r.get('status') == 'failed']
            failed_reports.extend(cur_failed_reports)
            save_fail_reports(cur_failed_reports, cur_step=wgd_status + '_all_solutions')
        save_fail_reports(failed_reports)

    if 'disambiguate' in which:
        logger.info('Starting KNN disambiguation of solutions with multiple paths')
        skip_existing = config['params'].get('skip_existing', False)
        full_paths_multiple_solutions_dirs=[os.path.join(results_dir, 'nowgd', 'full_paths_multiple_solutions'),
                                        os.path.join(results_dir, 'wgd', 'full_paths_multiple_solutions')]
        for wgd_status in ['nowgd', 'wgd']:
            if not os.path.exists(os.path.join(str(results_dir), wgd_status, 'full_paths_multiple_solutions')):
                logger.warning(f"Directory {os.path.join(str(results_dir), wgd_status, 'full_paths_multiple_solutions')} does not exist, skipping disambiguation for {wgd_status}")
                continue
            is_wgd = (wgd_status == 'wgd')
            cur_ids = [x.replace('.pickle', '')
                    for x in os.listdir(os.path.join(str(results_dir), wgd_status, 'full_paths_multiple_solutions'))]
            if selected_ids is not None:
                cur_ids = [x for x in cur_ids if x in selected_ids]
            def run_knn(cur_id):
                output_file = os.path.join(results_dir, wgd_status, 'knn_solved_chroms', f'{cur_id}.pickle')
                if skip_existing and os.path.exists(output_file):
                    logger.info(f"Skipping disambiguate for {cur_id} ({wgd_status}) as {output_file} exists.")
                    return {'status': 'skipped', 'cur_id': cur_id, 'step': 'disambiguate'}
                return solve_with_knn_wrapper(
                    output_file=output_file ,
                    cur_id=cur_id ,
                    is_wgd=is_wgd ,
                    chrom_segments_file=chrom_segments_file,
                    knn_train_data=None,
                    k=config['params']['knn_k'],
                    full_paths_multiple_solutions_dirs=full_paths_multiple_solutions_dirs,
                    save_all_scores=None,
                    perform_loh_checks=True,
                    single_width_bin=True
                )
            results = _run_batch(cur_ids, args.cores, f'Disambiguate solutions ({wgd_status})', run_knn, logger)
            cur_failed_reports = [r for r in results if isinstance(r, dict) and r.get('status') == 'failed']
            failed_reports.extend(cur_failed_reports)
            save_fail_reports(cur_failed_reports, cur_step=wgd_status + '_disambiguate')
        save_fail_reports(failed_reports)

    if 'large_chroms' in which:
        logger.info('Starting MCMC inference for large chromosomes with many events')
        skip_existing = config['params'].get('skip_existing', False)
        for wgd_status in ['nowgd', 'wgd']:
            if not os.path.exists(os.path.join(str(results_dir), wgd_status, 'chrom_data_large')):
                logger.warning(f"Directory {os.path.join(str(results_dir), wgd_status, 'chrom_data_large')} does not exist, skipping large chromosomes for {wgd_status}")
                continue
            is_wgd = (wgd_status == 'wgd')
            cur_ids = [x.replace('.pickle', '')
                    for x in os.listdir(os.path.join(str(results_dir), wgd_status, 'chrom_data_large'))]
            if selected_ids is not None:
                cur_ids = [x for x in cur_ids if x in selected_ids]

            def run_mcmc(cur_id):
                output_file = os.path.join(results_dir, wgd_status, 'mcmc_solved_chroms_large', f'{cur_id}.pickle')
                if skip_existing and os.path.exists(output_file):
                    logger.info(f"Skipping large_chroms for {cur_id} ({wgd_status}) as {output_file} exists.")
                    return {'status': 'skipped', 'cur_id': cur_id, 'step': 'large_chroms'}
                @timeout(config['params']['time_limit_mcmc'], mode="auto")
                def _solve_with_mcmc_wrapper(cur_id, skip_loh_check):
                    return solve_with_mcmc_wrapper(
                        output_file=output_file,
                        chrom_file=os.path.join(results_dir, wgd_status, 'chrom_data_large', f'{cur_id}.pickle'),
                        is_wgd=is_wgd,
                        chrom_segments_file=chrom_segments_file,
                        sv_data_file=config['params'].get('sv_data_file', None),
                        knn_train_data=None,
                        k=config['params']['knn_k'],
                        total_cn=total_cn,
                        save_all_scores=None,
                        n_iteration_scale=config['params']['mcmc_n_iterations_scale'],
                        log_progress=True,
                        fail_on_empty=False,
                        skip_loh_check=skip_loh_check
                    )
                try:
                    return _solve_with_mcmc_wrapper(cur_id, skip_loh_check=False)
                except FunctionTimeoutError as e:
                    logger.warning(f'MCMC solving for {cur_id} timed out after {config["params"]["time_limit_mcmc"]} seconds. Will rerun with LOH checks skipped. This might lead to inaccurate results! Consider increasing "time_limit_mcmc" in the config file.')
                    return _solve_with_mcmc_wrapper(cur_id, skip_loh_check=True)

            results = _run_batch(cur_ids, args.cores, f'Large chromosomes ({wgd_status})', run_mcmc, logger)
            cur_failed_reports = [r for r in results if isinstance(r, dict) and r.get('status') == 'failed']
            failed_reports.extend(cur_failed_reports)
            save_fail_reports(cur_failed_reports, cur_step=wgd_status + '_large_chroms')
        save_fail_reports(failed_reports)


    if 'combine' in which:
        logger.info('Starting combination of final events from all solving methods')
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
            output_dir=results_dir
        )

    save_fail_reports(failed_reports, logger=logger)
    logger.info(f'Done. Results are in {results_dir}')


def main_plotting(args):
    """Run plotting mode."""
    import pandas as pd
    
    # Load configuration
    spice.load_config(args.config_path)
    from spice import config
    from spice.data_loaders import load_final_events, resolve_data_file
    from spice.logging import configure_logging, get_logger
    from spice import plot as spice_plot
    from matplotlib import pyplot as plt

    if 'name' not in config or not config['name']:
        raise ValueError("Config file must specify a 'name' field.")
    
    # Create logger
    log_level = 'DEBUG' if args.debug else config['params'].get('logging_level', 'INFO')
    configure_logging(
        log_mode=args.log,
        log_dir=config['directories']['log_dir'],
        config_name=config['name'],
        level=log_level,
    )
    logger = get_logger('SPICE', spice_prefix=False)

    name = config['name']
    plots_base_dir = os.path.join(config['directories']['plot_dir'], name)
    if not os.path.exists(plots_base_dir):
        logger.info(f"Creating directory {plots_base_dir}")
        os.makedirs(plots_base_dir)

    logger.info('Running SPICE: Plotting Mode')
    logger.info(f'Plotting for project name {name} with config file {args.config_path}')

    # Validate plot target - this is now enforced by mutually_exclusive_group in argparse
    # Load required inputs
    chrom_segments_file = resolve_data_file()
    chrom_segments = pd.read_csv(
        chrom_segments_file, sep='\t', index_col=['sample_id', 'chrom', 'allele']).sort_index()
    final_events_df = load_final_events()

    if args.plot_sample is not None:
        cur_sample = args.plot_sample
        logger.info(f'Plotting inferred events for sample: {cur_sample}')
        fig, axs = spice_plot.plot_inferred_events_per_sample(
            cur_sample,
            chrom_segments,
            final_events_df,
            unit_size=args.plot_unit_size,
        )
        out_path = os.path.join(plots_base_dir, f'{cur_sample}_events{"_unit_size" if args.plot_unit_size else ""}.png')
        fig.savefig(out_path, bbox_inches='tight')
        logger.info(f'Saved plot to {out_path}')
    else:
        cur_id = args.plot_id
        logger.info(f'Plotting inferred events for id: {cur_id}')
        # Derive WIDTH_FULL from matplotlib defaults if not provided
        WIDTH_FULL = plt.rcParams.get('figure.figsize', (15, 5))[0]
        fig = spice_plot.plot_inferred_events_per_id(
            cur_id,
            chrom_segments,
            final_events_df,
            single_row=True,
            show_legend=True,
            figsize=(WIDTH_FULL, 1.25/5*WIDTH_FULL),
            lw=3,
            markersize=4
        )
        safe_id = cur_id.replace(':', '_')
        out_path = os.path.join(plots_base_dir, f'{safe_id}_events.png')
        fig.savefig(out_path, bbox_inches='tight')
        logger.info(f'Saved plot to {out_path}')

    logger.info('Done plotting.')


def main_loci_detection(args):
    """Run loci detection mode (de-novo)."""
    # Load configuration
    spice.load_config(args.config_path)
    from spice import config
    from spice.logging import configure_logging, get_logger

    if 'name' not in config or not config['name']:
        raise ValueError("Config file must specify a 'name' field.")
    
    # Create logger
    log_level = 'DEBUG' if args.debug else config['params'].get('logging_level', 'INFO')
    configure_logging(
        log_mode=args.log,
        log_dir=config['directories']['log_dir'],
        config_name=config['name'],
        level=log_level,
    )
    logger = get_logger('SPICE', spice_prefix=False)

    logger.info('Running SPICE: Loci Detection Mode (De-Novo)')
    logger.info(f'Project name: {config["name"]}')
    
    # Handle snakemake mode
    if args.snakemake:
        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
        snakefile = os.path.join(repo_root, 'Snakefile_loci_detection')
        if not os.path.exists(snakefile):
            raise FileNotFoundError(f"Snakefile_loci_detection not found at {snakefile}")

        cmd = [
            'snakemake',
            '-s', snakefile,
            '--rerun-triggers', 'mtime',
            '--verbose',
            '--configfile', args.config_path,
            '--config', f'config_path={args.config_path}',
        ]
        
        # Add number of jobs/cores
        if args.snakemake_mode == 'slurm':
            cmd.extend(['--slurm', '-j', str(args.snakemake_jobs), '--keep-going'])
        else:
            cmd.extend(['-c', str(args.snakemake_cores)])

        env = os.environ.copy()
        env['SPICE_CONFIG'] = os.path.abspath(args.config_path)
        logger.info(f'Running Snakemake loci detection workflow')
        subprocess.run(cmd, check=True, env=env)
        return
    
    logger.warning('Non-Snakemake mode for loci detection is not yet implemented.')
    logger.info('Please use --snakemake flag to run loci detection.')


def main():
    """Main CLI entry point for SPICE."""
    # Allow `spice --config <path> --snakemake` (default to event_inference mode)
    if '--snakemake' in sys.argv and not any(
        mode in sys.argv for mode in ['event_inference', 'plotting', 'loci_detection']
    ):
        sys.argv.insert(1, 'event_inference')

    parser = argparse.ArgumentParser(
        description='SPICE: Selection Patterns In somatic Copy-number Events',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Add --version flag at the top level
    parser.add_argument(
        '--version',
        action='version',
        version=f'%(prog)s {get_version()}'
    )
    
    parser.epilog = """
Examples:
  # Event inference
  spice event_inference --config <path/to/config>
  spice event_inference --config <path/to/config> --event-steps split all_solutions
  spice event_inference --config <path/to/config> --cores 8
  spice event_inference --config <path/to/config> --clean
  
  # Plotting
  spice plotting --config <path/to/config> --plot-sample "sample_1"
  spice plotting --config <path/to/config> --plot-id "sample_1:chr1:cn_a"
  
  # Loci detection (not yet implemented)
  spice loci_detection --config <path/to/config>
    """
    
    # Create subparsers for different modes
    subparsers = parser.add_subparsers(
        dest='mode',
        required=True,
        help='SPICE mode to run'
    )
    
    # Common arguments shared by all modes
    common_parser = argparse.ArgumentParser(add_help=False)
    common_parser.add_argument(
        '--config', '-c',
        required=True,
        type=str,
        dest='config_path',
        help='Path to a YAML config file to merge over defaults'
    )
    common_parser.add_argument(
        '--log',
        type=str,
        choices=['terminal', 'file', 'both'],
        default='terminal',
        help='Logging output mode: terminal (console only), file (log file only), or both (default: terminal)'
    )
    common_parser.add_argument(
        '--debug',
        action='store_true',
        help='Enable DEBUG logging globally, overriding config logging_level'
    )
    
    # ===== EVENT INFERENCE SUBPARSER =====
    parser_event = subparsers.add_parser(
        'event_inference',
        parents=[common_parser],
        help='Run event inference pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Infer discrete copy-number events from allele-specific profiles'
    )
    parser_event.add_argument(
        '--event-steps',
        nargs='+',
        default=argparse.SUPPRESS,
        help='Steps to run: preprocessing, split, all_solutions, disambiguate, large_chroms, combine (default: all). Use a trailing + (e.g., split+) to run that step and all subsequent steps.'
    )
    parser_event.add_argument(
        '--cores', '-j',
        type=int,
        default=None,
        help='Number of cores to use for parallel processing (default: 1)'
    )
    parser_event.add_argument(
        '--keep-old',
        action='store_true',
        help='Keep old intermediate files instead of overwriting them'
    )
    parser_event.add_argument(
        '--clean',
        action='store_true',
        help='Clean intermediate files and exit'
    )
    parser_event.add_argument(
        '--ids',
        type=str,
        default=None,
        help='Comma-separated list of sample IDs to process'
    )
    parser_event.add_argument(
        '--skip-preprocessing',
        action='store_true',
        help='Skip the extra preprocessing step that normally runs before split'
    )
    parser_event.add_argument(
        '--pre-unique-chroms',
        dest='pre_unique_chroms',
        action='store_true',
        help='Preprocessing: keep only unique chromosomes'
    )
    parser_event.add_argument(
        '--pre-skip-phasing',
        dest='pre_skip_phasing',
        action='store_true',
        help='Preprocessing: skip MEDICC2 phasing'
    )
    parser_event.add_argument(
        '--pre-skip-centromeres',
        dest='pre_skip_centromeres',
        action='store_true',
        help='Preprocessing: skip centromere binning'
    )
    parser_event.add_argument(
        '--snakemake',
        action='store_true',
        help='Run the event inference workflow using Snakemake instead of the Python runner'
    )
    parser_event.add_argument(
        '--snakemake-mode',
        choices=['local', 'slurm'],
        default='local',
        help='Snakemake execution mode: local or slurm (default: local)'
    )
    parser_event.add_argument(
        '--snakemake-jobs',
        type=int,
        default=250,
        help='Number of jobs for Snakemake on Slurm (-j, default: 250)'
    )
    parser_event.add_argument(
        '--snakemake-cores',
        type=int,
        default=1,
        help='Number of cores for local Snakemake execution (-c, default: 1)'
    )
    parser_event.add_argument(
        '--unlock',
        action='store_true',
        help='Unlock the Snakemake working directory and exit'
    )
    parser_event.set_defaults(func=main_event_inference)
    
    # ===== PLOTTING SUBPARSER =====
    parser_plot = subparsers.add_parser(
        'plotting',
        parents=[common_parser],
        help='Plot inferred events',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Generate visualizations of inferred copy-number events'
    )
    plot_group = parser_plot.add_mutually_exclusive_group(required=True)
    plot_group.add_argument(
        '--plot-sample',
        dest='plot_sample',
        type=str,
        help='Sample ID to plot'
    )
    plot_group.add_argument(
        '--plot-id',
        dest='plot_id',
        type=str,
        help='Chromosome allele ID to plot (format: sample:chr:cn_a|cn_b)'
    )
    parser_plot.add_argument(
        '--plot-unit-size',
        dest='plot_unit_size',
        action='store_true',
        help='Use unit_size for plotting events'
    )
    parser_plot.set_defaults(func=main_plotting)
    
    # ===== LOCI DETECTION SUBPARSER =====
    parser_loci = subparsers.add_parser(
        'loci_detection',
        parents=[common_parser],
        help='Detect recurrent copy-number loci (de-novo mode)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Identify recurrent copy-number loci across chromosomes using de-novo detection'
    )
    parser_loci.add_argument(
        '--snakemake',
        action='store_true',
        help='Run the loci detection workflow using Snakemake'
    )
    parser_loci.add_argument(
        '--snakemake-mode',
        choices=['local', 'slurm'],
        default='local',
        help='Snakemake execution mode: local or slurm (default: local)'
    )
    parser_loci.add_argument(
        '--snakemake-jobs',
        type=int,
        default=250,
        help='Number of jobs for Snakemake on Slurm (-j, default: 250)'
    )
    parser_loci.add_argument(
        '--snakemake-cores',
        type=int,
        default=1,
        help='Number of cores for local Snakemake execution (-c, default: 1)'
    )
    parser_loci.set_defaults(func=main_loci_detection)
    
    # Parse arguments and call the appropriate function
    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
