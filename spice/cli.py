#!/usr/bin/env python
"""Command-line interface for SPICE."""

import os
import shutil
import argparse

# Import base package only; defer submodule imports until after config is loaded
import spice

def main():
    """Main CLI entry point for SPICE."""
    # Defer logger creation until after config is loaded
    logger = None

    parser = argparse.ArgumentParser(
        description='SPICE: Selection Patterns In somatic Copy-number Events',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  spice --config <path/to/config>                               # Run all pipeline steps
  spice --config <path/to/config> split                         # Run only input splitting
  spice --config <path/to/config> all_solutions disambiguate    # Run path enumeration and kNN disambiguation
  spice --config <path/to/config> --cores 8                     # Run with 8 cores
  spice --config <path/to/config> plot --sample "sample_1"      # Run optional plotting for "sample_1"
    """
    )
    parser.add_argument(
        'steps',
        nargs='*',
        default=['all'],
        help='Steps to run: preprocessing, split, all_solutions, disambiguate, large_chroms, combine, plot (default: all). Use a trailing + (e.g., split+) to run that step and all subsequent steps.'
    )
    parser.add_argument(
        '--cores', '-j',
        required=False, 
        type=int, 
        default=None, 
        dest='cores',
        help='Number of cores to use for parallel processing (default: 1)'
    )
    parser.add_argument(
        '--config', '-c',
        required=True,
        type=str,
        default=None,
        dest='config_path',
        help='Path to a YAML config file to merge over defaults'
    )
    parser.add_argument(
        '--keep-old',
        action='store_true',
        help='Keep old intermediate files instead of overwriting them'
    )
    parser.add_argument(
        '--total_cn',
        action='store_true',
        help='Clean intermediate files'
    )
    parser.add_argument(
        '--clean',
        action='store_true',
        help='Clean intermediate files'
    )
    parser.add_argument(
        '--ids',
        type=str,
        default=None
    )
    parser.add_argument(
        '--log',
        type=str,
        choices=['terminal', 'file', 'both'],
        default='terminal',
        help='Logging output mode: terminal (console only), file (log file only), or both (default: terminal)'
    )
    parser.add_argument(
        '--debug',
        action='store_true',
        help='Enable DEBUG logging globally, overriding config logging_level'
    )
    parser.add_argument(
        '--skip-preprocessing',
        action='store_true',
        help='Skip the extra preprocessing step that normally runs before split'
    )
    # Preprocessing-specific arguments (optional)
    parser.add_argument(
        '--pre-unique-chroms',
        dest='pre_unique_chroms',
        action='store_true',
        help='Preprocessing: keep only unique chromosomes'
    )
    parser.add_argument(
        '--pre-skip-phasing',
        dest='pre_skip_phasing',
        action='store_true',
        help='Preprocessing: skip MEDICC2 phasing'
    )
    parser.add_argument(
        '--pre-skip-centromeres',
        dest='pre_skip_centromeres',
        action='store_true',
        help='Preprocessing: skip centromere binning'
    )
    # Plot-specific arguments
    parser.add_argument(
        '--plot-sample',
        dest='plot_sample',
        type=str,
        default=None,
        help='Sample ID to plot (for step: plot)'
    )
    parser.add_argument(
        '--plot-id',
        dest='plot_id',
        type=str,
        default=None,
        help='Chromosome allele ID to plot (format: sample:chr:cn_a|cn_b) (for step: plot)'
    )
    parser.add_argument(
        '--plot-unit-size',
        dest='plot_unit_size',
        action='store_true',
        help='Use unit_size for plotting events (applies to plot step)'
    )
    args = parser.parse_args()

    # Handle 'all' or empty arguments, and expand trailing + syntax (e.g., split+)
    step_order = ['preprocessing', 'split', 'all_solutions', 'disambiguate', 'large_chroms', 'combine']
    if args.steps == ['all'] or args.steps == []:
        which = step_order.copy()
    elif any('+' in x for x in args.steps):
        assert len(args.steps) == 1, 'Only one step with + is allowed'
        assert args.steps[0].endswith('+'), 'Only trailing + syntax is supported'
        which = step_order[step_order.index(args.steps[0][:-1]):]
    else:
        which = args.steps

    valid_steps = step_order + ['plot', 'all']
    invalid_steps = [step for step in which if step not in valid_steps]
    if invalid_steps:
        parser.error(f"Invalid step(s): {', '.join(invalid_steps)}. Valid steps are: preprocessing, split, all_solutions, disambiguate, large_chroms, combine, plot")

    # Load configuration before importing submodules that may read it
    spice.load_config(args.config_path)
    from spice import config
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
    logger.info(f'Running for project name {name} with config file {args.config_path}')

    if args.total_cn:
        raise NotImplementedError('--total-cn is not implemented yet')

    if args.clean:
        logger.info(f'Cleaning intermediate files at {results_dir}')
        for wgd in ['nowgd', 'wgd']:
            shutil.rmtree(os.path.join(results_dir, wgd), ignore_errors=True)
        return

    logger.info(f'Results will be stored in {results_dir}')
    logger.info(f'Running the following steps: {", ".join(which)}')

    selected_ids = args.ids.split(',') if args.ids is not None else None
    if selected_ids is not None:
        logger.info(f'Selecting only IDs: {selected_ids}')

    total_cn = config['input_files'].get('total_cn', False)
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
            total_cn=args.total_cn,
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
                    total_cn=total_cn,
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

    # Optional plotting step
    if 'plot' in which:
        import pandas as pd
        # Lazy import of plotting module to avoid heavy dependencies unless needed
        from spice import plot as spice_plot
        from matplotlib import pyplot as plt

        logger.info('Starting plotting of inferred events')

        # Validate plot target
        if not ((args.plot_sample is not None) ^ (args.plot_id is not None)):
            logger.error('Plot step requires exactly one of --plot-sample or --plot-id')
            return

        # Load required inputs
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

    save_fail_reports(failed_reports, logger=logger)

    logger.info(f'Done. Results are in {results_dir}')


if __name__ == '__main__':
    main()
