import pandas as pd
import numpy as np
from tqdm.auto import tqdm

from spice.utils import open_pickle, CALC_NEW
from spice.event_inference.SV import overlap_svs_with_events_df
from spice.event_inference.knn_graph import calc_event_distances


def calc_summary_from_events_df(events_df, chrom_segments):
    # chrom_segments is a DataFrame indexed by ['sample_id', 'chrom', 'allele'] with a 'cn' column
    # Build chrom_strings per chromosome/allele and include in summary
    chrom_strings = chrom_segments.groupby(['sample_id', 'chrom', 'allele'])['cn'].apply(lambda x: ''.join(x.astype(str)))
    chrom_strings = chrom_strings.to_frame('chrom_string')
    chrom_strings['id'] = chrom_strings.index.map(lambda x: ':'.join(x))
    chrom_strings = chrom_strings.reset_index().set_index('id')[['chrom_string']]

    summary_df = events_df.groupby('id')[['sample', 'chrom', 'allele', 'n_paths', 'events_per_chrom']].first()
    summary_df = summary_df.join(((chrom_segments.set_index('id')[['cn']] == 0).groupby('id')).sum().rename({'cn': 'lohs'}, axis=1), on='id')
    summary_df['sample_chrom_id'] = summary_df.index.map(lambda x: ':'.join(x.split(':')[:-1]))
    summary_df = summary_df.join((chrom_segments.groupby('id')['cn'].max() <= 1).to_frame('is_0_1'), on='id')
    summary_df = summary_df.join((summary_df.groupby('sample_chrom_id').size() == 2).to_frame('has_sister_id'), 
                on='sample_chrom_id')
    summary_df = summary_df.join(events_df.groupby('id')['wgd'].first().to_frame('WGD') != 'nowgd', on='id')
    summary_df = summary_df.join(chrom_strings, on='id')

    return summary_df


def overlap_final_events_with_svs(unique_events_df, sv_data, sv_matching_threshold):

    if sv_data is None or isinstance(sv_data, bool) or (isinstance(sv_data, str) and sv_data == ''):
        unique_events_df['SV_overlap'] = 0
    else:
        if isinstance(sv_data, str):
            sv_data = open_pickle(sv_data, fail_if_nonexisting=True)
        assert isinstance(sv_data, pd.DataFrame), type(sv_data)
        unique_events_df, sv_data = overlap_svs_with_events_df(
            unique_events_df, sv_data, relevant_chroms=None, threshold=sv_matching_threshold,
            filter_for_single_overlap=True, calc_new_verbose=False, verbose=False)

    return unique_events_df


def calculate_final_events_knn_score(unique_events_df, knn_train_data, knn_k, ignore_empty_train=False,
                                     clip_k=True):
    if knn_train_data is not None:
        if isinstance(knn_train_data, str):
            knn_train_data = open_pickle(knn_train_data, fail_if_nonexisting=True)
        elif isinstance(knn_train_data, dict):
            pass
        else:
            raise ValueError('knn_train_data must be a path to a pickle file or a dictionary')
        assert isinstance(knn_train_data, dict)
        event_distances = calc_event_distances(
            knn_train_data, unique_events_df, block_same_id=False,
            ks=[knn_k], show_progress=False, ignore_empty_train=ignore_empty_train, clip_k=clip_k)
        unique_events_df['knn_score'] = event_distances[:, 0]
    else:
        unique_events_df['knn_score'] = -1
    
    return unique_events_df


@CALC_NEW()
def calc_event_overlaps(final_events_df):

    assert (final_events_df.index == final_events_df.reset_index().index).all()
    final_events_df = final_events_df.copy()
    all_ids = final_events_df['id'].unique()

    final_events_df['nr_overlaps'] = 0
    final_events_df['nr_same_events'] = 0
    final_events_df['overlapping_events'] = None
    final_events_df['same_events'] = None
    for i in final_events_df.index:
        final_events_df.at[i, 'overlapping_events'] = []
        final_events_df.at[i, 'same_events'] = []

    for cur_id in tqdm(all_ids):
        cur_events = final_events_df.query('id == @cur_id')
        if len(cur_events) == 1:
            continue
        else:
            same_event = (
                (cur_events['start'].values == cur_events['start'].values[:, None]) &
                (cur_events['end'].values == cur_events['end'].values[:, None]))
            same_event[np.diag_indices_from(same_event)] = False
            cur_overlap = ~(
                (cur_events['end'].values < cur_events['start'].values[:, None]) |
                (cur_events['start'].values > cur_events['end'].values[:, None]))
            cur_overlap[np.diag_indices_from(cur_overlap)] = False
            cur_overlap = np.logical_and(cur_overlap, ~same_event)

            final_events_df.loc[cur_events.index, 'nr_overlaps'] = cur_overlap.sum(axis=0)
            for i, val in zip(cur_events.index, [list(np.where(x)[0]) for x in cur_overlap]):
                final_events_df.at[i, 'overlapping_events'] = cur_events.index[val]

            final_events_df.loc[cur_events.index, 'nr_same_events'] = same_event.sum(axis=0)
            for i, val in zip(cur_events.index, [list(np.where(x)[0]) for x in same_event]):
                final_events_df.at[i, 'same_events'] = cur_events.index[val]

    assert (final_events_df['overlapping_events'].map(len) == final_events_df['nr_overlaps']).all()
    assert (final_events_df['same_events'].map(len) == final_events_df['nr_same_events']).all()

    return final_events_df[['nr_same_events', 'same_events', 'overlapping_events', 'nr_overlaps']]