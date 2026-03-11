import os

from spice.data_loaders import format_chromosomes
import numpy as np
import pandas as pd
from tqdm import tqdm

from spice.utils import CALC_NEW, get_logger
from spice.data_loaders import load_chrom_lengths
from spice import config


logger = get_logger(__name__)
sv_matching_threshold = config['params']['sv_matching_threshold']


@CALC_NEW()
def overlap_svs_with_events_df(events_df, sv_data, relevant_chroms=None, threshold=sv_matching_threshold,
                               verbose=True, filter_for_single_overlap=False):

    events_df = events_df.copy()
    sv_data = sv_data.copy()
    events_df['SV_overlap'] = 0
    sv_data['event_overlap'] = 0

    if relevant_chroms is None:
        relevant_chroms = np.intersect1d(events_df['chrom_id'].unique(),
        # relevant_chroms = np.intersect1d(events_df.query('n_paths > 1')['chrom_id'].unique(),
                                         sv_data['chrom_id'].unique())
        # relevant_chroms = events_df['chrom_id'].unique()
    logger.debug(f"overlap svs with events for {len(relevant_chroms)} chroms")

    for cur_chrom_id in tqdm(relevant_chroms, disable=not verbose):
        cur_sv_data = sv_data.query('chrom_id == @cur_chrom_id and (svclass == "DUP" or svclass == "DEL")').copy()
        cur_events_df = events_df.query('chrom_id == @cur_chrom_id').copy()
        
        cur_events_df, cur_sv_data = overlap_svs_with_events_df_single(
            cur_events_df, cur_sv_data, sv_matching_threshold=threshold,
            filter_for_single_overlap=filter_for_single_overlap)
        events_df.loc[cur_events_df.index, 'SV_overlap'] = cur_events_df['SV_overlap']
        sv_data.loc[cur_sv_data.index, 'event_overlap'] = cur_sv_data['event_overlap']

    return events_df, sv_data


def overlap_svs_with_events_df_single(cur_events_df, cur_sv_data, sv_matching_threshold=sv_matching_threshold,
                                      filter_for_single_overlap=False):

    cur_events_df['SV_overlap'] = 0
    cur_sv_data['event_overlap'] = 0

    # check that all values of events_df.index are unique
    assert len(cur_events_df.index) == len(cur_events_df.index.unique()), 'not all index values are unique'

    # Note: if it fails here, it is probably because the index of events_df is not unique -> reset_index()
    events_sv_overlap_gain = np.logical_and(
        np.abs(cur_sv_data.query('svclass == "DUP"')['start'].values - cur_events_df.query('type == "gain"')['start'].values[:, np.newaxis]) < sv_matching_threshold,
        np.abs(cur_sv_data.query('svclass == "DUP"')['end'].values - cur_events_df.query('type == "gain"')['end'].values[:, np.newaxis]) < sv_matching_threshold)
    if 0 not in events_sv_overlap_gain.shape and np.sum(events_sv_overlap_gain) > 0:
        cur_sv_data.loc[cur_sv_data.query('svclass == "DUP"').index, 'event_overlap'] += np.sum(events_sv_overlap_gain, axis=0)
        if filter_for_single_overlap:
            cur_mask = np.sum(events_sv_overlap_gain, axis=0) == 1
            if np.any(cur_mask):
                events_sv_overlap_gain = events_sv_overlap_gain[:, cur_mask]
                cur_events_df.loc[cur_events_df.query('type == "gain"').index, 'SV_overlap'] += np.sum(events_sv_overlap_gain, axis=1)
        else:
            cur_events_df.loc[cur_events_df.query('type == "gain"').index, 'SV_overlap'] += np.sum(events_sv_overlap_gain, axis=1)

    events_sv_overlap_loss = np.logical_and(
        np.abs(cur_sv_data.query('svclass == "DEL"')['start'].values - cur_events_df.query('type == "loss"')['start'].values[:, np.newaxis]) < sv_matching_threshold,
        np.abs(cur_sv_data.query('svclass == "DEL"')['end'].values - cur_events_df.query('type == "loss"')['end'].values[:, np.newaxis]) < sv_matching_threshold)
    if 0 not in events_sv_overlap_loss.shape and np.sum(events_sv_overlap_loss) > 0:
        cur_sv_data.loc[cur_sv_data.query('svclass == "DEL"').index, 'event_overlap'] += np.sum(events_sv_overlap_loss, axis=0)
        if filter_for_single_overlap:
            cur_mask = np.sum(events_sv_overlap_loss, axis=0) == 1
            if np.any(cur_mask):
                events_sv_overlap_loss = events_sv_overlap_loss[:, cur_mask]
                cur_events_df.loc[cur_events_df.query('type == "loss"').index, 'SV_overlap'] += np.sum(events_sv_overlap_loss, axis=1)
        else:
            cur_events_df.loc[cur_events_df.query('type == "loss"').index, 'SV_overlap'] += np.sum(events_sv_overlap_loss, axis=1)

    return cur_events_df, cur_sv_data
