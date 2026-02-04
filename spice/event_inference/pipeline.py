import os
from collections import namedtuple
import itertools

import pandas as pd
import numpy as np
import fire

from spice import config
from spice.utils import (
    save_pickle, open_pickle, chrom_id_from_id,
    calc_telomere_bound_whole_arm_whole_chrom)
from spice.logging import get_logger, log_debug
from spice.event_inference.events_from_graph import (
    full_paths_from_graph_with_sv, create_events_df_from_single_path_solution)
from spice.event_inference.knn_graph import solve_with_knn
from spice.event_inference.mcmc_for_large_chroms import mcmc_event_selection, create_best_events_df_from_mcmc
from spice.event_inference.data_structures import ChromData
from spice.pipeline_postprocessing import calc_summary_from_events_df
from spice.event_analysis.final_events import classify_event_position


logger = get_logger('spice.pipeline')


EVENT_DF_COLUMNS = [
    'sample', 'chrom', 'id', 'width', 'start', 'end', 'type', 'telomere_bound', 'whole_chrom', 'whole_arm',
    'diff', 'wgd', 'allele', 'chrom_id', 'events_per_chrom', 'chrom_length', 'solved', 'n_paths']
EVENT_DF_DTYPES = {
    'sample': "str",
    'chrom': "str",
    'id': "str",
    'width': "int64",
    'start': "int64",
    'end': "int64",
    'type': "str",
    'telomere_bound': "bool",
    'whole_chrom': "bool",
    'whole_arm': "bool",
    'diff': "str",
    'wgd': "str",
    'allele': "str",
    'chrom_id': "str",
    'events_per_chrom': "int64",
    'chrom_length': "int64",
    'solved': "str",
    'n_paths': "int64",
}


def full_paths_from_graph_with_sv_wrapper(
        cur_id, is_wgd, sv_data_file, chrom_segments_file, chrom_file,
        sv_matching_threshold=10, without_sv_output_dir=False,
        all_loh_solutions=True, total_cn=False,
        skip_loh_checks=False, use_cache=True, save_output=True):
    
    # Implemente here to debug the creation of the fail report
    # if cur_id == "RPelvicLNMet_A12D-0020_CRUK_PC_0020_M3_DEBUG:chr1:cn_a":
    #     raise ValueError("Debug ID")
    
    log_debug(logger, f"Full path solving for {cur_id} ({'WGD' if is_wgd else 'noWGD'})")
    chrom_id = chrom_id_from_id(cur_id)
    if sv_data_file is None or sv_data_file == '' or isinstance(sv_data_file, bool):
        cur_sv_data = None
        log_debug(logger, 'No SV data provided, so skipping SV overlap')
    else:
        log_debug(logger, f'Loading SV data from {sv_data_file}')
        sv_data = open_pickle(sv_data_file)
        assert isinstance(sv_data, pd.DataFrame)
        cur_sv_data = sv_data.query('chrom_id == @chrom_id')
    chrom_segments = pd.read_csv(chrom_segments_file, sep='\t', index_col=['sample_id', 'chrom', 'allele'])
    chrom = open_pickle(chrom_file)

    full_paths = full_paths_from_graph_with_sv(
        cur_id, is_wgd, cur_sv_data, chrom_segments, chrom, path_limit=None,
        sv_matching_threshold=sv_matching_threshold, total_cn=total_cn,
        all_loh_solutions=all_loh_solutions, use_cache=use_cache, skip_loh_checks=skip_loh_checks)

    if len(full_paths.events) == 0 or len(full_paths.solutions) == 0:
        raise ValueError(f'No events (n = {len(full_paths.events)}) or no solutions (n = {len(full_paths.solutions)}) found')

    if full_paths.solved in ['sv', 'unamb']:
        full_paths_dir = 'full_paths_without_sv_single_solution' if without_sv_output_dir else 'full_paths_single_solution'
        log_debug(logger, f"Full path is completely solved with {full_paths.solved}")
        assert len(full_paths.solutions) == 1, f"solved_chrom.solutions: {full_paths.solutions}"
        output = create_events_df_from_single_path_solution(
            full_paths, cur_id, chrom_segments=chrom_segments, create_full=True, calc_telomere_bound=False)
    else:
        full_paths_dir = 'full_paths_without_sv_multiple_solutions' if without_sv_output_dir else 'full_paths_multiple_solutions'
        output = full_paths


    if save_output:
        log_debug(logger, f'saving file to {config["directories"]["results_dir"]}/{config["name"]}/{"wgd" if is_wgd else "nowgd"}/{full_paths_dir}/{cur_id}.pickle')
        save_pickle(output, os.path.join(config['directories']['results_dir'], config["name"], 'wgd' if is_wgd else 'nowgd', full_paths_dir, f'{cur_id}.pickle'))
    else:
        return output


def create_knn_train_data(output_file, knn_train_data_ext=None, full_paths_single_solution_dirs=None,
                          knn_train_dataset='sv_unamb'):
    if knn_train_data_ext is not None and knn_train_data_ext != '':
        log_debug(logger, f'external knn_train_data provided at {knn_train_data_ext}, so using it')
        knn_train_data = open_pickle(knn_train_data_ext, fail_if_nonexisting=True)
    else:
        log_debug(logger, f'no external knn_train_data provided, so creating new one from files from {full_paths_single_solution_dirs}')
        assert full_paths_single_solution_dirs is not None
        if isinstance(full_paths_single_solution_dirs, str):
            full_paths_single_solution_dirs = full_paths_single_solution_dirs.split(' ')

        unamb_events = []
        for cur_dir in full_paths_single_solution_dirs:
            if not os.path.exists(cur_dir):
                logger.warning(f'{cur_dir} does not exist, skipping')
                continue
            for cur_file in os.listdir(cur_dir):
                if not cur_file.endswith('.pickle'):
                    continue
                cur_events = open_pickle(os.path.join(cur_dir, cur_file), fail_if_nonexisting=True)
                unamb_events.append(cur_events)
            
        unamb_events = pd.concat(unamb_events).reset_index(drop=True)
        assert unamb_events['solved'].isin(['sv', 'unamb']).all(), 'All events must be either unambiguous or solved by SVs'
        if knn_train_dataset not in ['sv', 'unamb', 'sv_unamb', 'unamb_sv']:
            raise ValueError(f"Invalid knn_train_dataset: {knn_train_dataset}. Must be in ['sv', 'unamb', 'sv_unamb', 'unamb_sv']")
        if knn_train_dataset == 'sv':
            log_debug(logger, 'selecting only SV solved events')
            unamb_events = unamb_events.query('solved == "sv"').copy()
        elif knn_train_dataset == 'unamb':
            log_debug(logger, 'selecting only unambiguous events')
            unamb_events = unamb_events.query('solved == "unamb"').copy()
        else:
            log_debug(logger, 'selecting all events (SV solved and unambiguous)')

        # recalculate this here because of bug in the pipeline (should be unnecessary as of 20.06.2024)
        unamb_events[['telomere_bound', 'whole_arm', 'whole_chrom']] = np.stack(
            calc_telomere_bound_whole_arm_whole_chrom(unamb_events), axis=1)

        log_debug(logger, f"Found {unamb_events['id'].nunique()} unique IDs with {len(unamb_events)} events")
        pre_post_wgds = ['nowgd', 'pre', 'post', 'all']

        knn_train_data = dict()
        for cur_type, telomere_bound, pre_post_wgd in itertools.product(['gain', 'loss'], [True, False], pre_post_wgds):        
            if pre_post_wgd == 'all':
                cur_mask_train = np.logical_and.reduce([
                    (unamb_events['type'] == cur_type),
                    (unamb_events['telomere_bound'] == telomere_bound),
                    (~unamb_events['whole_chrom']),
                    (~unamb_events['whole_arm'])
                ])
            else:
                cur_mask_train = np.logical_and.reduce([
                    (unamb_events['type'] == cur_type),
                    (unamb_events['telomere_bound'] == telomere_bound),
                    (unamb_events['wgd'] == pre_post_wgd),
                    (~unamb_events['whole_chrom']),
                    (~unamb_events['whole_arm'])
                ])
            knn_train_data[(cur_type, telomere_bound, pre_post_wgd)] = (
                unamb_events.loc[cur_mask_train, 'width'].values,
                unamb_events.loc[cur_mask_train, 'chrom_length'].values,
                unamb_events.loc[cur_mask_train, 'id'].values,
                )
    
    save_pickle(knn_train_data, output_file)


def solve_with_knn_wrapper(cur_id, chrom_segments_file, full_paths_multiple_solutions_dirs, is_wgd, knn_train_data=None,
                           k=250, single_width_bin=True, perform_loh_checks=False,
                           save_all_scores=None, output_file=None):
    
    if isinstance(full_paths_multiple_solutions_dirs, str):
        full_paths_multiple_solutions_dirs = full_paths_multiple_solutions_dirs.split(' ')

    full_paths = open_pickle(os.path.join(full_paths_multiple_solutions_dirs[1 if is_wgd else 0], f'{cur_id}.pickle'), fail_if_nonexisting=True)
    assert getattr(full_paths, '__class__', None).__name__ == 'FullPaths', type(full_paths)
    # assert isinstance(full_paths, FullPaths) 

    cur_chrom_segments = pd.read_csv(chrom_segments_file, sep='\t', index_col=['sample_id', 'chrom', 'allele']).query('id == @cur_id')

    log_debug(logger, f"KNN-graph solving for {cur_id} ({'WGD' if is_wgd else 'noWGD'}). CN profile: {full_paths.cn_profile}, n_solutions: {full_paths.n_solutions}, n_events: {full_paths.n_events}")

    if full_paths.solved in ['sv', 'unamb']:
        raise ValueError('full_paths already solved')
    if len(full_paths.events) == 0:
        raise ValueError('full_paths does not have any events')

    if knn_train_data is None:
        knn_train_data = open_pickle(config['input_files']['knn_train'], fail_if_nonexisting=True)
    elif not isinstance(knn_train_data, dict):
        log_debug(logger, f'Loading knn_train_data from {knn_train_data}')
        knn_train_data = open_pickle(knn_train_data, fail_if_nonexisting=True)
        assert isinstance(knn_train_data, dict)
    final_events_df = solve_with_knn(
        full_paths,
        k=k,
        cur_chrom_segments=cur_chrom_segments,
        knn_train_data=knn_train_data,
        perform_loh_checks=perform_loh_checks,
        wgd_split=True,
        log10_distances=True,
        save_all_scores=save_all_scores,
        single_width_bin=single_width_bin)

    final_events_df = final_events_df[EVENT_DF_COLUMNS]

    if output_file is not None:
        save_pickle(final_events_df, output_file)
    else:
        return final_events_df


def solve_with_mcmc_wrapper(
        chrom_file, chrom_segments_file, sv_data_file, is_wgd, knn_train_data=None, k=250, wgd_split=True,
        output_file=None, sv_matching_threshold=10, n_iterations=None, n_iteration_scale=100, perform_loh_checks=False,
        min_T=1, max_T=-6, swap_event_based_on_score=True, check_all_loh_solutions=False, total_cn=False,
        verbose=False, save_all_scores=None, log_progress=False, show_progress=False, fail_on_empty=True,
        skip_loh_check=False):
    assert (n_iterations is not None) ^ (n_iteration_scale is not None), 'Either n_iterations or n_iteration_scale must be provided'

    chrom_data = open_pickle(chrom_file, fail_if_nonexisting=True)
    assert isinstance(chrom_data, ChromData)
    cur_id = chrom_data.id

    log_debug(logger, f"MCMC solving for {cur_id} ({'WGD' if is_wgd else 'noWGD'})")
    assert is_wgd == chrom_data.has_wgd, f'is_wgd ({is_wgd}) does not match chrom_data.has_wgd ({chrom_data.has_wgd})'
    cur_chrom_segments = pd.read_csv(chrom_segments_file, sep='\t', index_col=['sample_id', 'chrom', 'allele']).query('id == @cur_id')

    if n_iterations is None:
        n_iterations = int(chrom_data.n_events * n_iteration_scale)

    if knn_train_data is None:
        knn_train_data =open_pickle(config['input_files']['knn_train'], fail_if_nonexisting=True)
    elif not isinstance(knn_train_data, dict):
        log_debug(logger, f'Loading knn_train_data from {knn_train_data}')
        knn_train_data = open_pickle(knn_train_data, fail_if_nonexisting=True)
        assert isinstance(knn_train_data, dict)

    if sv_data_file is None or sv_data_file == '' or isinstance(sv_data_file, bool):
        cur_sv_data = None
        log_debug(logger, 'No SV data provided, so skipping SV overlap')
    else:
        log_debug(logger, f'Loading SV data from {sv_data_file}')
        chrom_id = chrom_id_from_id(cur_id)
        sv_data = open_pickle(sv_data_file, fail_if_nonexisting=True)
        assert isinstance(sv_data, pd.DataFrame)
        cur_sv_data = sv_data.query('chrom_id == @chrom_id')
   
    mcmc_result = mcmc_event_selection(
                        cur_id=cur_id,
                        chrom_data=chrom_data,
                        chrom_segments=cur_chrom_segments,
                        sv_data=cur_sv_data,
                        acceptance_temp=None,
                        n_iterations=n_iterations,
                        knn_graph_k=k,
                        sv_matching_threshold=sv_matching_threshold,
                        knn_train_data=knn_train_data,
                        show_progress=show_progress,
                        simulated_annealing=True,
                        min_T=min_T,
                        max_T=max_T,
                        wgd_split=wgd_split,
                        loh_time_limit=None,
                        swap_event_based_on_score=swap_event_based_on_score,
                        perform_loh_checks=perform_loh_checks,
                        check_all_loh_solutions=check_all_loh_solutions,
                        total_cn=total_cn,
                        log_progress=log_progress,
                        calc_new_filename=None,
                        calc_new_verbose=False,
                        skip_loh_check=skip_loh_check
                        )
    
    # if it fails with SVs, rerun without them
    if mcmc_result is None and (cur_sv_data is not None and len(cur_sv_data) > 0):
        cur_sv_data = None
        log_debug(logger, 'No solution found with SVs, rerunning without SVs')
        mcmc_result = mcmc_event_selection(
                        cur_id=cur_id,
                        chrom_data=chrom_data,
                        chrom_segments=cur_chrom_segments,
                        sv_data=cur_sv_data,
                        acceptance_temp=None,
                        n_iterations=n_iterations,
                        knn_graph_k=k,
                        sv_matching_threshold=sv_matching_threshold,
                        knn_train_data=knn_train_data,
                        show_progress=show_progress,
                        simulated_annealing=True,
                        min_T=min_T,
                        max_T=max_T,
                        wgd_split=wgd_split,
                        loh_time_limit=None,
                        swap_event_based_on_score=swap_event_based_on_score,
                        check_all_loh_solutions=check_all_loh_solutions,
                        total_cn=total_cn,
                        perform_loh_checks=perform_loh_checks,
                        log_progress=log_progress,
                        calc_new_filename=None,
                        calc_new_verbose=False)

    if mcmc_result is not None:
        final_events_df = create_best_events_df_from_mcmc(mcmc_result, cur_chrom_segments, chrom_data, has_wgd=chrom_data.has_wgd)  
        final_events_df = final_events_df[EVENT_DF_COLUMNS]
    else:
        if fail_on_empty:
            raise ValueError('MCMC run failed, no solution was found.')
        final_events_df = None
   
    if save_all_scores is not None:
        log_debug(logger, f'Saving all scores to {save_all_scores}')
        save_pickle(mcmc_result, save_all_scores)

    if output_file is not None:
        log_debug(logger, f'Saving final events to {output_file}')
        save_pickle(final_events_df, output_file)

    return final_events_df

def combine_final_events(solved_dirs, chrom_segments_file=None, sv_data=None,
                         sv_matching_threshold=10, knn_train_data=None, knn_k=250, output_dir=None):

    if isinstance(solved_dirs, str):
        solved_dirs = solved_dirs.split(' ')

    chrom_segments = pd.read_csv(chrom_segments_file, sep='\t', index_col=['sample_id', 'chrom', 'allele'])

    solved_events = []
    for cur_dir_i, cur_dir in enumerate(solved_dirs):
        if not os.path.exists(cur_dir):
            log_debug(logger, f'{cur_dir} does not exist, skipping')
            continue
        cur_dir_name = '/'.join(cur_dir.split('/')[-2:])
        log_debug(logger, f'Directory {cur_dir_i+1}/{len(solved_dirs)}: {cur_dir}')
        n_files = len(os.listdir(os.path.join(cur_dir)))
        for i, cur_file in enumerate(os.listdir(os.path.join(cur_dir))):
            if not cur_file.endswith('.pickle'):
                continue
            log_debug(logger, f'Directory {cur_dir_i+1}/{len(solved_dirs)}, file {i+1}/{n_files} ({100*(i+1)/n_files:.2f}%): {cur_file} from {cur_dir_name}')
            cur_events = open_pickle(os.path.join(cur_dir, cur_file), fail_if_nonexisting=True)
            if cur_events is None:
                log_debug(logger, f'Data is None, skipping')
                continue
            solved_events.append(cur_events[EVENT_DF_COLUMNS].values)
    log_debug(logger, f'Found {len(solved_events)} solved IDs')

    if len(solved_events) == 0:
        raise ValueError('No solved events found in the provided directories')
    final_events_df = pd.DataFrame(np.concatenate(solved_events, axis=0), columns=EVENT_DF_COLUMNS)
    log_debug(logger, f'Found a total of {len(final_events_df)} events in the final df')

    # missing_telomere_bound has to be before types are assigned
    missing_telomere_bound = np.logical_or(
        (final_events_df[['telomere_bound', 'whole_arm', 'whole_chrom']] == None).any(axis=1),
        (final_events_df[['telomere_bound', 'whole_arm', 'whole_chrom']].isna()).any(axis=1)
    )
    final_events_df.loc[missing_telomere_bound, ['telomere_bound', 'whole_arm', 'whole_chrom']] = (
        np.stack(calc_telomere_bound_whole_arm_whole_chrom(final_events_df.loc[missing_telomere_bound]), axis=1))

    final_events_df = final_events_df[EVENT_DF_COLUMNS]
    final_events_df = final_events_df.astype(dtype=EVENT_DF_DTYPES)

    final_events_df = final_events_df.join((final_events_df.groupby('id')['wgd'].first() != 'nowgd').to_frame('has_wgd'), on='id')
    final_events_df['loh_shortened_event'] = final_events_df['diff'].map(lambda x: '0' in x[x.find('1'):x.rfind('1')])
    final_events_df['pos'] = classify_event_position(final_events_df)

    # if knn_train_data is None:
    #     open_pickle(config['input_files']['knn_train'], fail_if_nonexisting=True)
    # unique_events_df = final_events_df.drop_duplicates(['id', 'diff', 'wgd']).copy().reset_index(drop=True)
    # log_debug(logger, 'Overlapping events with SVs')
    # unique_events_df = overlap_final_events_with_svs(unique_events_df, sv_data, sv_matching_threshold)
    # log_debug(logger, 'Calculating final event knn scores')
    # unique_events_df = calculate_final_events_knn_score(unique_events_df, knn_train_data, knn_k)
    # final_events_df = final_events_df.join(unique_events_df.set_index(['id', 'diff', 'wgd'])[['SV_overlap', 'knn_score']],
    #                                        on=['id', 'diff', 'wgd']).reset_index(drop=True)
    # log_debug(logger, f'Joined SV overlaps and knn scores into main df, new length = {len(final_events_df)}')
    # final_events_df = final_events_df.join(
    #     (chrom_segments.groupby('id')['end'].max() - chrom_segments.groupby('id')['start'].min()).to_frame('actual_chrom_length'),
    #     on='id')
    # final_events_df = (final_events_df
    #     .join((final_events_df.set_index('id')['wgd'] == 'pre').groupby('id').sum().to_frame('n_pre'), on='chrom')
    #     .join((final_events_df.set_index('id')['wgd'] == 'post').groupby('id').sum().to_frame('n_post'), on='chrom')
    #     )

    log_debug(logger, f'After postprocessing, found a total of {len(final_events_df)} events in the final df')
    log_debug(logger, 'Calculating summary')
    summary_df = calc_summary_from_events_df(final_events_df, chrom_segments)

    log_debug(logger, f'Saving final files to {output_dir}')
    if output_dir is not None:
        final_events_df.to_csv(os.path.join(output_dir, 'final_events.tsv'), sep='\t', index=False)
        summary_df.to_csv(os.path.join(output_dir, 'summary.tsv'), sep='\t', index=True)
        # save_pickle(final_events_df, os.path.join(output_dir, 'final_events.pickle'))
        # save_pickle(summary_df, os.path.join(output_dir, 'summary.pickle'))

    return final_events_df, summary_df


def final_quality_control(events_file, log_dir=None):

    events_df = pd.read_csv(events_file, sep='\t', dtype={'id': str, 'diff': str})
    events_df['loh_shortened_event'] = events_df['diff'].map(lambda x: '0' in x[x.find('1'):x.rfind('1')])

    file_overview = []
    for folder in ['chrom_data_full', 'chrom_data_large', 'full_paths_created', 'full_paths_multiple_solutions',
                   'full_paths_single_solution', 'knn_solved_chroms', 'mcmc_solved_chroms_large', 'mcmc_solved_chroms_full']:
        cur_files = []
        for wgd in ['nowgd', 'wgd']:
            cur_dir = os.path.join(os.path.dirname(events_file), wgd, folder)
            if not os.path.exists(cur_dir):
                cur_files.append(pd.DataFrame(columns=[folder]))
                continue
            cur_ids = [x.replace('.pickle', '').replace('.placeholder', '') for x in os.listdir(cur_dir) if x.endswith('.pickle') or x.endswith('.placeholder')]
            cur_files.append(pd.DataFrame(True, index=cur_ids, columns=[folder]))
        if len(cur_files) == 0:
            cur_files.append(pd.DataFrame(False, columns=[folder]))
            continue 
        cur = pd.concat(cur_files, axis=1)
        cur = ~cur.isna()
        file_overview.append(cur.any(axis=1).to_frame(folder))
    file_overview = pd.concat(file_overview, axis=1)
    # fillna does not work with bools
    file_overview = ~file_overview.isna()
    file_overview.to_csv(os.path.join(os.path.dirname(events_file), 'file_overview.tsv'), sep='\t', index=True)

    Test = namedtuple('Test', ['test', 'type', 'message'])
    tests = [
        Test(not events_df.isna().any().any(),
             'final_df', f'some columns have NAs: {events_df.columns[events_df.isna().any()].values}'),
        Test((events_df.query('solved=="single"')['events_per_chrom'] == 1).all(),
             'final_df', 'not all single have single event'),
        Test(len(events_df) == events_df.drop_duplicates('id')['events_per_chrom'].sum(),
            'final_df', "length of df ({len(events_df)}) does not match expected nr of events ({events_df.drop_duplicates('id')['events_per_chrom'].sum()}) "
            f" for IDs {(events_df.groupby('id').size() != events_df.groupby('id')['events_per_chrom'].first()).to_frame('test').query('test').index.values}"),
        Test(len(events_df.index) == len(events_df.index.unique()),
             'final_df', 'not all index values are unique'),
        Test(events_df.query('solved=="sv"').groupby('id')['SV_overlap'].any().all(),
            'final_df', 'not all sv solved events have sv overlap'),
        Test(np.isfinite(events_df['knn_score'].max()),
             'final_df', 'knn score is not finite'),
        Test(np.isin(events_df['type'], ['gain', 'loss']).all(),
             'final_df', 'not all events are gain or loss'),
        Test(events_df.loc[events_df['diff'].str.startswith('1') & events_df['diff'].str.endswith('1'), 'whole_chrom'].all(),
             'final_df', 'incorrect whole-chrom assignment'),
        Test(events_df.loc[events_df['diff'].str.startswith('1') | events_df['diff'].str.endswith('1'), 'telomere_bound'].all(),
             'final_df', 'incorrect telomere-bound assignment'),
        Test(events_df.eval('loh_shortened_event or (end - start) == width').all(),
             'final_df', 'not all events have end - start == width'),
        Test(events_df.loc[events_df['diff'].str.endswith('1')].eval('end == chrom_length').all(),
             'final_df', 'not all ..1 events end at chromosome length'),
        Test(events_df.eval('whole_chrom or ((width != actual_chrom_length) or (width != chrom_length))').all(),
             'final_df', 'non-whole-chrom events span the whole chromosome'),
        Test(events_df.loc[events_df['diff'].str.startswith('1')].eval('start == 0 or ((chrom == "chr13") or (chrom == "chr14") or (chrom == "chr15") or (chrom == "chr21") or (chrom != "chr22"))').all(),
                'final_df', 'not all 1.. events start at 0'),
        Test(events_df.loc[events_df['diff'].str.count('1') == events_df['diff'].str.len()].eval('width == actual_chrom_length').all(),
                'final_df', 'not all whole chrom events have width == chrom_length'),
        Test(events_df.eval('actual_chrom_length <= chrom_length').all(),
                'final_df', 'not all actual_chrom_length <= chrom_length'),
        Test(events_df.eval('not whole_chrom or (whole_chrom and telomere_bound)').all(),
                'final_df', 'not all whole_chrom events are telomere bound'),
        Test(events_df.eval('not whole_arm or (whole_arm and telomere_bound)').all(),
                'final_df', 'not all whole_arm events are telomere bound'),
        Test(events_df.eval('not (whole_chrom and whole_arm)').all(),
                'final_df', 'some events are both whole_chrom and whole_arm'),
        
        Test(np.isin(events_df['id'], file_overview.index.values).all(),
            'final_df', f'not all ids from final_df are in file_overview: {np.setdiff1d(events_df["id"], file_overview.index.values)}'),
        # The below errouneously fails for MCMC that time out
        # Test(np.isin(file_overview.index.values, events_df['id']).all(),
        #     'file_overview', f'not all ids from file_overview are in final_df: {np.setdiff1d(file_overview.index.values, events_df["id"])}'),
    ]

    # tuples of length 2: if file exists in any of the dirs in the first entry it hast to exist in all of
    # the first entries and cannot exist in any of the dirs of the second entry, and vice verse
    dir_groups = [
        (['chrom_data_full'], 
         ['chrom_data_large', 'mcmc_solved_chroms_large']),
        (['full_paths_multiple_solutions', 'knn_solved_chroms'], 
         ['chrom_data_large', 'mcmc_solved_chroms_large']),
        (['full_paths_multiple_solutions', 'knn_solved_chroms'],
        ['full_paths_single_solution'])
        ]
    
    for dir_group in dir_groups:
        tests.extend([
            Test(file_overview.loc[file_overview[dir_group[0]].any(axis=1)][dir_group[0]].all(axis=1).all(),
                 'file_overview', f'Some files are missing for: {dir_group[0]}'),
            Test(file_overview.loc[file_overview[dir_group[1]].any(axis=1)][dir_group[1]].all(axis=1).all(),
                 'file_overview', f'Some files are missing for: {dir_group[1]}'),
            Test(not file_overview.loc[file_overview[dir_group[0]].any(axis=1)][dir_group[1]].any(axis=1).any(),
                 'file_overview', f'Some files are in both groups for: {dir_group}'),
            Test(not file_overview.loc[file_overview[dir_group[1]].any(axis=1)][dir_group[0]].any(axis=1).any(),
                 'file_overview', f'Some files are in both groups for: {dir_group}')
        ])

    # as full path creation can fail, have to check that files from chrom_data_full are solved either by knn or mcmc
    tests.append(Test(file_overview.loc[file_overview['chrom_data_full'], ['full_paths_created', 'mcmc_solved_chroms_full']].any(axis=1).all(),
                        'file_overview', 'Some files are missing for chrom_data_full'))


    passed = True
    for test in tests:
        if not test.test:
            logger.error(f'{test.type}: {test.message}')
            passed = False
    
    assert passed, 'Quality control failed'
    log_debug(logger, f"passed all {len(tests)} tests")

    if log_dir is not None:
        import subprocess

        command = f"grep 'WARNING' {log_dir}/*"
        try:
            result = subprocess.run(command, shell=True, check=True, text=True, capture_output=True)
            warnings = [w.replace(log_dir+'/', ' ') for w in result.stdout.split('\n') if len(w) > 0]
            log_debug(logger, f"Found {len(warnings)} WARNINGS in the logs ({log_dir})")
            with open(os.path.join(os.path.dirname(events_file), 'log_warnings.txt'), 'w') as f:
                f.write('\n'.join(warnings))
        except subprocess.CalledProcessError as e:
            logger.warning(f"Error while checking log directory: {e}")

        
if __name__ == '__main__':
  fire.Fire()
