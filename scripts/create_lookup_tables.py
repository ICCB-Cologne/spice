"""
Create lookup tables for single and multiple solution full paths.
These are made up of common copy number profiles found in chrom_strings.csv. 
During inference, these lookup tables can be used to quickly retrieve precomputed full paths
instead of recalculating them.

Usage:
    python scripts/create_lookup_tables.py <chrom_strings_path> [--output-dir OUTPUT_DIR]

Example:
    python scripts/create_lookup_tables.py results/analysis/lookup_events/chrom_strings.csv
    python scripts/create_lookup_tables.py data/chrom_strings.csv --output-dir objects/
"""

import os
import logging
import argparse
import numpy as np
import pandas as pd
from collections import Counter
from tqdm import tqdm

from spice.event_inference.data_structures import Diff, FullPaths
from spice.event_inference.events_from_graph import (
    get_events_from_graph_wgd,
    get_events_diff_from_coords_wgd,
    get_starts_and_ends,
    _full_paths_implementation_nowgd,
    loh_filters_for_graph_result_diffs,
    loh_filters_for_graph_result_diffs_wgd
)
from spice import config
from spice.utils import log_debug, save_pickle, get_logger

# Set up logging
logger = get_logger('create_lookup_tables', load_config=True)
logging.getLogger('spice.events_from_graph').setLevel(logging.WARNING)


def create_full_path_from_profile(cn_profile, is_wgd):
    """Create FullPaths object from CN profile and WGD status."""
    if is_wgd:
        # Generate candidate paths and diffs
        paths = get_events_from_graph_wgd(cn_profile)
        diffs = get_events_diff_from_coords_wgd(
            paths, cn_profile, lexsort_diffs=True, filter_missed_lohs=False
        )
        # Apply LOH filters (matching full_paths_from_graph_with_sv behavior)
        if 0 in cn_profile:
            diffs = loh_filters_for_graph_result_diffs_wgd(
                diffs,
                cn_profile,
                single_time_limit=config['params'].get('time_limit_loh_filters', 60),
                return_all_solutions=False,
                total_cn=False,
                shuffle_diffs=True,
            )
        if len(diffs) == 0:
            raise NotImplementedError('No viable solutions (diffs) found after LOH filtering')
        # Convert to Diff objects
        diffs = [
            [Diff(diff=''.join(map(str, np.abs(x))), is_gain=x.max() == 1, wgd="pre") for x in diff[0]]
            + [Diff(diff=''.join(map(str, np.abs(x))), is_gain=x.max() == 1, wgd="post") for x in diff[1]]
            for diff in diffs
        ]
    else:
        starts, ends = get_starts_and_ends(cn_profile, prior_profile=None, loh_adjust=True)
        diffs = _full_paths_implementation_nowgd(None, starts, ends, starts, cn_profile)
        # Apply LOH filters for non-WGD if LOH present
        if 0 in cn_profile:
            diffs = loh_filters_for_graph_result_diffs(
                diffs,
                cn_profile,
                single_time_limit=config['params'].get('time_limit_loh_filters', 60),
                total_cn=False,
                return_all_solutions=False,
                shuffle_diffs=True,
            )
        if len(diffs) == 0:
            raise NotImplementedError('No viable solutions (diffs) found after LOH filtering')
        diffs = [
            [Diff(diff=''.join(map(str, np.abs(x))), is_gain=x.max() == 1, wgd="nowgd") for x in diff]
            for diff in diffs
        ]
    
    # Create unique events dictionary
    unique_events = {i: d for i, d in enumerate(set(item for sublist in diffs for item in sublist))}
    unique_events_reversed = {v: k for k, v in unique_events.items()}
    diffs = [[unique_events_reversed[event] for event in diff] for diff in diffs]
    
    # Create solutions and remove duplicates
    solutions = [Counter(diff) for diff in diffs]
    unique_solutions = [Counter({k: v for k, v in x}) for x in {frozenset(c.items()) for c in solutions}]
    
    full_path = FullPaths(
        id='test:chr1:cn_a', 
        sample='test', 
        chrom='chr1', 
        allele='cn_a',
        cn_profile=cn_profile,
        solutions=unique_solutions, 
        events=unique_events,
        n_solutions=len(unique_solutions), 
        n_events=len(diffs[0]),
        is_wgd=is_wgd, 
        solved=None
    )
    
    return full_path


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Create lookup tables for single and multiple solution full paths.'
    )
    parser.add_argument(
        'chrom_strings_path',
        type=str,
        help='Path to the chrom_strings.csv file'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default=None,
        help='Output directory for lookup tables (default: objects/ in repo root)'
    )
    parser.add_argument(
        '--n-single',
        type=int,
        default=None,
        help='Number of single solution entries to include (default: all)'
    )
    parser.add_argument(
        '--n-multiple',
        type=int,
        default=500,
        help='Number of multiple solution entries to include (default: 500)'
    )
    parser.add_argument(
        '--debug',
        action='store_true',
        help='Enable DEBUG logging globally, overriding config logging_level'
    )
    args = parser.parse_args()
    
    if args.debug:
        logger.setLevel(logging.DEBUG)

    # Validate input file exists
    chrom_strings_path = args.chrom_strings_path
    assert os.path.exists(chrom_strings_path), f"chrom_strings.csv not found at {chrom_strings_path}"
    
    # Determine output directory
    if args.output_dir is None:
        script_dir = os.path.dirname(__file__)
        lookup_dir = os.path.join(script_dir, '..', 'objects')
    else:
        lookup_dir = args.output_dir
    os.makedirs(lookup_dir, exist_ok=True)
    logger.info(f"Using output directory: {lookup_dir}")
    
    # Load chrom_strings
    chrom_strings = pd.read_csv(chrom_strings_path, sep='\t', dtype={'cn': str})
    logger.info(f"Loaded chrom_strings from {chrom_strings_path}")
    chrom_strings['cn'] = chrom_strings['cn'].astype(str)
    
    # Validate input
    required_columns = ['solved', 'cn', 'WGD', 'count']
    for col in required_columns:
        assert col in chrom_strings.columns, f"Missing required column: {col}, found columns: {chrom_strings.columns.tolist()}"
    
    # Create single solution lookup table
    logger.info("Creating single solution lookup table...")
    single_solution_full_paths = dict()
    
    single_solutions_data = chrom_strings.query('solved == "unamb"').reset_index(drop=True)
    if args.n_single is not None:
        single_solutions_data = single_solutions_data.iloc[:args.n_single]
    for i, row in tqdm(single_solutions_data.iterrows(), total=len(single_solutions_data), desc="Single solutions", disable=args.debug):
        cn_profile = np.fromiter(row['cn'], dtype=int)
        log_debug(logger, f"Single solution {i}/{len(single_solutions_data)} ({100*i/len(single_solutions_data):.2f}%): {cn_profile}, {'WGD' if row['WGD'] else 'noWGD'}")
        full_path = create_full_path_from_profile(cn_profile, row['WGD'])
        single_solution_full_paths[(row['WGD'], row['cn'])] = full_path
    
    # Create multiple solutions lookup table
    logger.info("Creating multiple solutions lookup table...")
    multiple_solutions_full_paths = dict()
    
    multiple_solutions_data = chrom_strings.query('solved == "knn"').sort_values('count', ascending=False).reset_index(drop=True)
    if args.n_multiple is not None:
        multiple_solutions_data = multiple_solutions_data.iloc[:args.n_multiple]
    for i, row in tqdm(multiple_solutions_data.iterrows(), total=len(multiple_solutions_data), desc="Multiple solutions", disable=args.debug):
        cn_profile = np.fromiter(row['cn'], dtype=int)
        log_debug(logger, f"Multiple solutions {i}/{len(multiple_solutions_data)} ({100*i/len(multiple_solutions_data):.2f}%): {cn_profile}, {'WGD' if row['WGD'] else 'noWGD'}")
        full_path = create_full_path_from_profile(cn_profile, row['WGD'])
        multiple_solutions_full_paths[(row['WGD'], row['cn'])] = full_path
    
    # Save lookup tables    
    single_path = os.path.join(lookup_dir, 'lookup_table_single_solution_full_paths.pickle')
    save_pickle(single_solution_full_paths, single_path)
    logger.info(f"Saved single solution lookup table to {single_path}")
    
    multiple_path = os.path.join(lookup_dir, 'lookup_table_multiple_solutions_full_paths.pickle')
    save_pickle(multiple_solutions_full_paths, multiple_path)
    logger.info(f"Saved multiple solutions lookup table to {multiple_path}")
    
    logger.info(f"Created {len(single_solution_full_paths)} single solution entries")
    logger.info(f"Created {len(multiple_solutions_full_paths)} multiple solution entries")


if __name__ == '__main__':
    main()
