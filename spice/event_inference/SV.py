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
def load_pcawg_sv_data(sample_labels, chrom_offset=None, sv_size_filter=None, filter_dup_del=False,
                       sv_data_dir='/Users/tom/phd/projects/medicc/data/PCAWG/SV/all'):
    
    if chrom_offset is None:
        chrom_offset = load_chrom_lengths().to_frame()
        chrom_offset['chrom_offset'] = np.append(0, chrom_offset.cumsum().values[:-1])
        chrom_offset = chrom_offset.drop(['chrX', 'chrY'])[['chrom_offset']]

    sv_data = []
    for sample in tqdm(sample_labels):
        if os.path.exists(os.path.join(sv_data_dir, '{}.pcawg_consensus_1.6.161116.somatic.sv.bedpe.gz'.format(sample))):
            cur_sv_data = pd.read_csv(os.path.join(sv_data_dir, '{}.pcawg_consensus_1.6.161116.somatic.sv.bedpe.gz'.format(sample)), 
                    compression='infer', sep='\t')
            cur_sv_data['sample'] = sample
            cur_sv_data = cur_sv_data.loc[cur_sv_data['chrom1']==cur_sv_data['chrom2']]
            cur_sv_data = cur_sv_data.loc[~cur_sv_data['chrom1'].isin(['X', 'Y'])]
            if len(cur_sv_data) == 0:
                continue

            if sv_size_filter is not None:
                cur_sv_data = cur_sv_data.loc[cur_sv_data['end2']-cur_sv_data['start1']>sv_size_filter]

            cur_sv_data['chrom'] = format_chromosomes(cur_sv_data['chrom1'])
            cur_sv_data = cur_sv_data.set_index(['chrom', 'start1', 'end2'])
            cur_sv_data = cur_sv_data.join(chrom_offset, how='inner')
            cur_sv_data['start_pos'] = cur_sv_data.eval('start1+chrom_offset')
            cur_sv_data['end_pos'] = cur_sv_data.eval('end2+chrom_offset')
            if filter_dup_del:
                cur_sv_data = cur_sv_data.loc[cur_sv_data['svclass'].isin(['DUP', 'DEL'])]
            
            if len(cur_sv_data) == 0:
                continue

            sv_data.append(cur_sv_data)
    sv_data = pd.concat(sv_data).reset_index()
    sv_data['width'] = sv_data.eval('end2-start1')
    sv_data['chrom_id'] = sv_data['sample'] + ':' + sv_data['chrom'].astype(str)

    return sv_data


@CALC_NEW()
def load_jabba_pcawg_sv_data(sample_labels, chrom_offset=None, sv_size_filter=None, filter_dup_del=False,
                             sv_data_dir='/Users/tom/phd/projects/signatures/data/jabba_pcawg/del_dup_calls',
                             _pcawg_sample_sheet='/Users/tom/phd/projects/cn_data/data/PCAWG/PCAWG_Sample_Sheet.tsv'):
    
    if not filter_dup_del:
        raise NotImplementedError('filter_dup_del=False is not implemented and does not make sense')

    _pcawg_sample_sheet = pd.read_csv(_pcawg_sample_sheet, sep='\t')
    _pcawg_id_dict = _pcawg_sample_sheet.set_index('icgc_donor_id')['aliquot_id'].to_dict()
    

    if chrom_offset is None:
        chrom_offset = load_chrom_lengths().to_frame()
        chrom_offset['chrom_offset'] = np.append(0, chrom_offset.cumsum().values[:-1])
        chrom_offset = chrom_offset.drop(['chrX', 'chrY'])['chrom_offset']

    sv_data = []
    for sample in tqdm(sample_labels):
        if os.path.exists(os.path.join(sv_data_dir, f'{sample}.txt')) and sample in _pcawg_id_dict:
            cur_sv_data = pd.read_csv(os.path.join(sv_data_dir, f'{sample}.txt'), sep='\t')
            cur_sv_data['sample'] = _pcawg_id_dict[sample]
            cur_sv_data = cur_sv_data.loc[cur_sv_data['chrom1']==cur_sv_data['chrom2']]
            cur_sv_data = cur_sv_data.loc[~cur_sv_data['chrom1'].isin(['X', 'Y'])]
            if len(cur_sv_data) == 0:
                continue

            cur_sv_data['svclass'] = cur_sv_data['class'].str.split('-').str[0]
            cur_sv_data = cur_sv_data.loc[cur_sv_data['svclass'].isin(['DUP', 'DEL'])]
            cur_sv_data.loc[cur_sv_data['svclass']=='DUP', ['start1', 'end2']] = cur_sv_data.loc[cur_sv_data['svclass']=='DUP', ['start2', 'end1']].values

            if sv_size_filter is not None:
                cur_sv_data = cur_sv_data.loc[cur_sv_data['end2']-cur_sv_data['start1']>sv_size_filter]

            cur_sv_data['chrom'] = format_chromosomes(cur_sv_data['chrom1'])
            cur_sv_data = cur_sv_data.set_index(['chrom', 'start1', 'end2'])
            cur_sv_data = cur_sv_data.join(chrom_offset, how='inner')
            cur_sv_data['start_pos'] = cur_sv_data.eval('start1+chrom_offset')
            cur_sv_data['end_pos'] = cur_sv_data.eval('end2+chrom_offset')
            
            if len(cur_sv_data) == 0:
                continue

            sv_data.append(cur_sv_data)
    sv_data = pd.concat(sv_data).reset_index()
    sv_data['width'] = sv_data.eval('end2-start1')
    sv_data['chrom_id'] = sv_data['sample'] + ':' + sv_data['chrom'].astype(str)

    return sv_data


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
        np.abs(cur_sv_data.query('svclass == "DUP"')['start1'].values - cur_events_df.query('type == "gain"')['start'].values[:, np.newaxis]) < sv_matching_threshold,
        np.abs(cur_sv_data.query('svclass == "DUP"')['end2'].values - cur_events_df.query('type == "gain"')['end'].values[:, np.newaxis]) < sv_matching_threshold)
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
        np.abs(cur_sv_data.query('svclass == "DEL"')['start1'].values - cur_events_df.query('type == "loss"')['start'].values[:, np.newaxis]) < sv_matching_threshold,
        np.abs(cur_sv_data.query('svclass == "DEL"')['end2'].values - cur_events_df.query('type == "loss"')['end'].values[:, np.newaxis]) < sv_matching_threshold)
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


@CALC_NEW()
def get_sv_supported_chains(events_df):

    selected_events_df = events_df.query('n_paths > 1 and type != "loh" and SV_overlap != -1').copy()

    selected_events_df['supported_chain'] = False
    for cur_id in tqdm(selected_events_df.query('SV_overlap > 0')['id'].unique()):

        cur_supported_diffs = selected_events_df.query('id == @cur_id and SV_overlap > 0')['diff'].unique()
        supported_chains = (selected_events_df
                            .query('id == @cur_id and SV_overlap > 0')
                            .drop_duplicates(['chain_nr', 'diff'])
                            .groupby('chain_nr')
                            .apply(lambda x: len(np.setdiff1d(cur_supported_diffs, x['diff'])))
                            .to_frame('missing_supported_diffs')
                            .query('missing_supported_diffs == 0')
                            .index.values)
        selected_events_df.loc[selected_events_df.query('id == @cur_id and chain_nr in @supported_chains').index, 'supported_chain'] = True

    selected_events_df = selected_events_df.drop('has_supported_chain', axis=1, errors='ignore').join(
        selected_events_df.groupby('id')['supported_chain'].any().to_frame('has_supported_chain'),
            on='id', how='left')

    return selected_events_df



@CALC_NEW()
def overlap_svs_with_events_df_legacy(events_df, sv_data, relevant_chroms=None, threshold=sv_matching_threshold, verbose=True):
    raise DeprecationWarning("use overlap_svs_with_events_df instead")

    events_df = events_df.copy()
    sv_data = sv_data.copy()
    events_df['SV_overlap'] = 0
    sv_data['event_overlap'] = 0

    if relevant_chroms is None:
        relevant_chroms = np.intersect1d(events_df.query('n_paths > 1')['chrom_id'].unique(),
                                         sv_data['chrom_id'].unique())
        # relevant_chroms = events_df['chrom_id'].unique()
    logger.debug(f"overlap svs with events for {len(relevant_chroms)} chroms")

    for cur_chrom_id in tqdm(relevant_chroms, disable=not verbose):
        cur_sv_data = sv_data.query('chrom_id == @cur_chrom_id and (svclass == "DUP" or svclass == "DEL")')
        cur_events_df = events_df.query('chrom_id == @cur_chrom_id')

        # Note: if it fails here, it is probably because the index of events_df is not unique -> reset_index()
        events_sv_overlap_gain = np.logical_and(
            np.abs(cur_sv_data.query('svclass == "DUP"')['start1'].values - cur_events_df.query('type == "gain"')['start'].values[:, np.newaxis]) < threshold,
            np.abs(cur_sv_data.query('svclass == "DUP"')['end2'].values - cur_events_df.query('type == "gain"')['end'].values[:, np.newaxis]) < threshold)
        if 0 not in events_sv_overlap_gain.shape and np.sum(events_sv_overlap_gain) > 0:
            sv_data.loc[cur_sv_data.query('svclass == "DUP"').index, 'event_overlap'] += np.sum(events_sv_overlap_gain, axis=0)
            events_df.loc[cur_events_df.query('type == "gain"').index, 'SV_overlap'] += np.sum(events_sv_overlap_gain, axis=1)

        events_sv_overlap_loss = np.logical_and(
            np.abs(cur_sv_data.query('svclass == "DEL"')['start1'].values - cur_events_df.query('type == "loss"')['start'].values[:, np.newaxis]) < threshold,
            np.abs(cur_sv_data.query('svclass == "DEL"')['end2'].values - cur_events_df.query('type == "loss"')['end'].values[:, np.newaxis]) < threshold)
        if 0 not in events_sv_overlap_loss.shape and np.sum(events_sv_overlap_loss) > 0:
            sv_data.loc[cur_sv_data.query('svclass == "DEL"').index, 'event_overlap'] += np.sum(events_sv_overlap_loss, axis=0)
            events_df.loc[cur_events_df.query('type == "loss"').index, 'SV_overlap'] += np.sum(events_sv_overlap_loss, axis=1)
        

    return events_df, sv_data

